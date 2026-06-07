# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Reproducibly add XAGE2 to the bundled CTA table and recompute the
parameterized inclusion columns.

This is a one-off curation script (kept in-repo for auditability) that:

1. Builds the XAGE2 (``ENSG00000155622``) evidence row from HPA
   ``rna_tissue_consensus.tsv`` using the same in-repo enrichment functions
   (:func:`tsarina.tiers.enrich_rna_per_tissue`, :func:`tsarina.tiers.assign_all_axes`)
   that generate the rest of the table.  XAGE2 is a CTpedia (CT12.2) gene that
   was missing from the source universe (see tsarina#79).
2. Recomputes ``passes_filters`` for every row from the documented adaptive
   reproductive-fraction rule, using the (now 0.97) Missing/Uncertain threshold
   in :data:`tsarina.tissues.HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS`.  This flips
   CT83 (KK-LC-1) and PRM3 to pass.
3. Recomputes ``never_expressed`` from the parameterized
   :data:`tsarina.tissues.HPA_EXPRESSION_FLOOR_NTPM` floor.

Run ``python scripts/add_xage2.py --check`` to validate the reimplemented
``passes_filters`` / ``never_expressed`` rules reproduce the shipped columns at
the *old* 0.98 threshold (proves the rule is faithful) and preview the XAGE2
row + flipped genes, without writing.  Run without ``--check`` to write the
updated CSV.

HPA version note: existing rows' per-tissue/fraction columns reflect the HPA
consensus release used when the table was first built; XAGE2's row is computed
from the current release.  The deflated reproductive fraction for placenta-
dominant genes is robust to this drift (validated against XAGE3 == 0.9922,
XAGE5 == 1.0).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from cta_sources import rna_tissue_consensus_path  # noqa: E402

from tsarina import tiers  # noqa: E402
from tsarina.tissues import (  # noqa: E402
    CORE_REPRODUCTIVE_TISSUES,
    HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS,
    HPA_EXPRESSION_FLOOR_NTPM,
)

CSV_PATH = REPO_ROOT / "tsarina" / "data" / "cancer-testis-antigens.csv"

XAGE2_ENSG = "ENSG00000155622"
XAGE2_IDENTITY = {
    "Symbol": "XAGE2",
    "Aliases": "XAGE2B",
    "Full_Name": None,
    "Function": "cancer-testis antigen",
    "Ensembl_Gene_ID": XAGE2_ENSG,
    "source_databases": "CTpedia",
    "Canonical_Transcript_ID": "ENST00000286049",
    "biotype": "protein_coding",
    # No HPA IHC protein data ("Not detected").
    "protein_reproductive": "no data",
    "protein_thymus": "no data",
    "protein_reliability": "no data",
    "protein_strict_expression": "no data",
    # No public MHC-ligand MS evidence yet (recomputed live from hitlist at
    # panel time); matches the default carried by other no-MS rows.
    "ms_restriction": "NO_MS_DATA",
}

_NO_PROTEIN = {"no data", "nan", ""}


def _has_protein(reliability) -> bool:
    return str(reliability).strip().lower() not in _NO_PROTEIN


def _protein_reproductive(value) -> bool:
    return str(value).strip().lower() == "true"


def passes_filters_rule(row: pd.Series, missing_threshold: float) -> bool:
    """Documented adaptive reproductive-fraction filter.

    Reproduces the bundled ``passes_filters`` column exactly at the original
    0.98 threshold (validated over all rows).
    """
    if str(row["biotype"]) != "protein_coding":
        return False
    reliability = str(row["protein_reliability"]).strip()
    has_protein = _has_protein(reliability)
    if has_protein and not _protein_reproductive(row["protein_reproductive"]):
        # Protein detected in a non-reproductive tissue -> always fails.
        return False
    tier = reliability if has_protein else "Missing"
    threshold = HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS.get(tier, missing_threshold)
    if tier in ("Uncertain", "Missing"):
        threshold = missing_threshold
    frac = row["rna_deflated_reproductive_frac"]
    if pd.isna(frac):
        return False
    return bool(float(frac) >= threshold)


def never_expressed_rule(row: pd.Series, floor: float) -> bool:
    """``never_expressed`` == no protein data AND max RNA nTPM < floor."""
    if _has_protein(row["protein_reliability"]):
        return False
    max_ntpm = row["rna_max_ntpm"]
    if pd.isna(max_ntpm):
        return True
    return bool(float(max_ntpm) < floor)


def _fraction(ntpm: dict[str, float], allowed: frozenset[str], deflate: bool) -> float:
    vals = {t: (max(0.0, v - 1.0) if deflate else v) for t, v in ntpm.items()}
    total = sum(vals.values())
    if total <= 0:
        return 1.0 if deflate else 0.0
    repro = sum(v for t, v in vals.items() if t in allowed)
    return repro / total


def build_xage2_row(consensus_path: Path, template_columns: list[str]) -> dict:
    """Construct the full XAGE2 evidence row from HPA consensus data."""
    consensus = pd.read_csv(consensus_path, sep="\t")
    sub = consensus[consensus["Gene"] == XAGE2_ENSG]
    if sub.empty:
        raise SystemExit(f"XAGE2 ({XAGE2_ENSG}) not found in {consensus_path}")
    ntpm = {t.strip().lower(): float(v) for t, v in zip(sub["Tissue"], sub["nTPM"])}

    core = CORE_REPRODUCTIVE_TISSUES
    core_thymus = frozenset(core | {"thymus"})

    row = dict.fromkeys(template_columns)
    row.update(XAGE2_IDENTITY)

    # RNA fraction columns (CORE reproductive set; validated against XAGE3).
    row["rna_reproductive_frac"] = round(_fraction(ntpm, core, deflate=False), 4)
    row["rna_reproductive_and_thymus_frac"] = round(_fraction(ntpm, core_thymus, deflate=False), 4)
    deflated = _fraction(ntpm, core, deflate=True)
    row["rna_deflated_reproductive_frac"] = round(deflated, 4)
    row["rna_deflated_reproductive_and_thymus_frac"] = round(
        _fraction(ntpm, core_thymus, deflate=True), 4
    )

    # Percentile filter booleans use the deflated reproductive fraction.
    for pct, thr in [("80", 0.80), ("90", 0.90), ("95", 0.95), ("97", 0.97), ("98", 0.98), ("99", 0.99)]:
        row[f"rna_{pct}_pct_filter"] = bool(deflated >= thr)

    row["rna_max_ntpm"] = round(max(ntpm.values()), 1)
    row["rna_thymus"] = bool(ntpm.get("thymus", 0.0) >= 1.0)

    # Per-tissue + safety columns, restriction axes -- via the in-repo generators.
    seed = pd.DataFrame([row])
    seed = tiers.enrich_rna_per_tissue(seed, str(consensus_path))
    # rna_reproductive: reproductive-restricted at RNA (no somatic tissue >= 1 nTPM).
    seed["rna_reproductive"] = seed["rna_somatic_detected_count"].fillna(0).astype(int).eq(0)
    seed = tiers.assign_all_axes(seed)

    out = seed.iloc[0].to_dict()
    out["passes_filters"] = passes_filters_rule(
        seed.iloc[0], HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS["Missing"]
    )
    out["never_expressed"] = never_expressed_rule(seed.iloc[0], HPA_EXPRESSION_FLOOR_NTPM)
    return {c: out.get(c) for c in template_columns}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--consensus",
        default=None,
        help="Local HPA rna_tissue_consensus path; default downloads + caches from HPA.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate rules reproduce shipped columns at 0.98; preview only, do not write.",
    )
    args = parser.parse_args()

    df = pd.read_csv(CSV_PATH)
    columns = list(df.columns)

    # ── Validate the reimplemented rules reproduce the shipped columns ──
    repro_pf = df.apply(lambda r: passes_filters_rule(r, 0.98), axis=1)
    shipped_pf = df["passes_filters"].astype(bool)
    pf_mismatch = int((repro_pf != shipped_pf).sum())

    repro_ne = df.apply(lambda r: never_expressed_rule(r, HPA_EXPRESSION_FLOOR_NTPM), axis=1)
    shipped_ne = df["never_expressed"].astype(str).str.lower() == "true"
    ne_mismatch = int((repro_ne != shipped_ne).sum())

    print(f"passes_filters reproduction @0.98: {len(df) - pf_mismatch}/{len(df)} match")
    print(f"never_expressed reproduction @floor={HPA_EXPRESSION_FLOOR_NTPM}: "
          f"{len(df) - ne_mismatch}/{len(df)} match")
    if pf_mismatch or ne_mismatch:
        raise SystemExit("Rule reproduction failed; aborting (would corrupt the table).")

    # ── Recompute at the new threshold + add XAGE2 ──
    new_threshold = HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS["Missing"]
    new_pf = df.apply(lambda r: passes_filters_rule(r, new_threshold), axis=1)
    flipped = df.loc[new_pf != shipped_pf, "Symbol"].tolist()
    print(f"\nThreshold {new_threshold}: passes_filters flips fail->pass: {flipped}")

    xage2 = build_xage2_row(rna_tissue_consensus_path(args.consensus), columns)
    print("\nXAGE2 row:")
    for key in ["Symbol", "Ensembl_Gene_ID", "source_databases", "rna_placenta_ntpm",
                "rna_max_somatic_tissue", "rna_max_somatic_ntpm", "rna_deflated_reproductive_frac",
                "passes_filters", "never_expressed", "restriction", "restriction_confidence",
                "safety_flags"]:
        print(f"  {key}: {xage2.get(key)}")

    if args.check:
        print("\n--check: no files written.")
        return

    df["passes_filters"] = new_pf
    df["never_expressed"] = repro_ne
    df = pd.concat([df, pd.DataFrame([xage2])[columns]], ignore_index=True)
    df.to_csv(CSV_PATH, index=False)
    print(f"\nWrote {len(df)} rows to {CSV_PATH}")


if __name__ == "__main__":
    main()
