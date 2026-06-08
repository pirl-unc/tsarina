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

"""Regenerate every RNA-derived column in the bundled CTA table from current HPA.

The bundled ``cancer-testis-antigens.csv`` accumulated **mixed HPA provenance**:
the filter-driving ``rna_deflated_reproductive_frac`` is (mostly) current, but
the per-tissue detail columns (``rna_max_ntpm``, ``rna_max_somatic_*``,
``rna_*_max_ntpm``) are from an older HPA release.  Evidence: 27 rows reference
brain subregions (pons/thalamus/medulla oblongata/white matter) absent from the
current 51-tissue consensus, and ``rna_max_ntpm`` is below a per-tissue value in
253 rows (impossible within one snapshot).  See the tissue-derivation audit.

This script re-derives **all RNA and protein/IHC columns** for **every existing
row** from a single pinned HPA release (``tsarina reference``; default v23 --
the newest release serving both ``rna_tissue_consensus`` and ``normal_tissue``),
using the same generators the per-gene ``add_cta_gene.py`` uses, so the whole
table is internally consistent and reproducible from one release.  Identity,
Aliases, source_databases, and MS columns are **preserved** verbatim (MS comes
from hitlist, not HPA).

It is **safe by default**: it writes a side-by-side ``*.regen.csv`` and prints a
membership/column delta report, but does NOT touch the shipped CSV unless
``--apply`` is passed.

Usage::

    python scripts/regenerate_table.py            # report + write sidecar, no overwrite
    python scripts/regenerate_table.py --apply     # overwrite the bundled CSV
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from add_cta_gene import (  # noqa: E402
    _fraction,
    never_expressed_rule,
    passes_filters_rule,
)
from cta_sources import normal_tissue_path, rna_tissue_consensus_path  # noqa: E402

from tsarina import tiers  # noqa: E402
from tsarina.tiers import _ALL_REPRODUCTIVE  # noqa: E402
from tsarina.tissues import (  # noqa: E402
    CORE_REPRODUCTIVE_TISSUES,
    HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS,
    HPA_EXPRESSION_FLOOR_NTPM,
)

#: HPA IHC levels that count as "protein detected" in a tissue.
_PROTEIN_DETECTED_LEVELS = frozenset({"Low", "Medium", "High"})
#: Reproductive (+thymus) tissues, lowercased, for the protein restriction call.
_PROTEIN_REPRODUCTIVE = frozenset(t.lower() for t in _ALL_REPRODUCTIVE)

#: Genes whose v23 HPA IHC is treated as **unreliable** (forced to "no data"),
#: so they fall back to the RNA-only restriction call.  These are sequence-near-
#: identical CT/paralog antigens whose shared antibody cross-reacts: the somatic
#: protein "detected" by HPA is at Low level in tissues where the gene's RNA is
#: ~0 nTPM (heart/glandular cells), the hallmark of cross-reactivity.  Curated
#: override -- keyed by (unversioned) Ensembl gene ID.  See the regeneration
#: discussion in the PR; revisit if HPA ships gene-specific antibodies.
_CROSS_REACTIVE_IHC: frozenset[str] = frozenset(
    {
        "ENSG00000176746",  # MAGEB6  -- kidney IHC, RNA 0; Enhanced but RNA-discordant
        "ENSG00000155622",  # XAGE2   -- scattered Low glandular IHC, RNA ~0 (placenta-restricted)
        "ENSG00000278085",  # CT45A8  -- heart Low IHC, RNA 0; CT45A1 paralog (shared antibody)
        "ENSG00000270946",  # CT45A9  -- heart Low IHC, RNA 0; CT45A1 paralog (shared antibody)
    }
)

CSV_PATH = REPO_ROOT / "tsarina" / "data" / "cancer-testis-antigens.csv"

#: RNA columns this script re-derives.  Everything else (identity, protein/IHC,
#: MS, Aliases, source_databases, ...) is preserved from the existing row.
_RNA_FRACTION_COLUMNS = [
    "rna_reproductive_frac",
    "rna_reproductive_and_thymus_frac",
    "rna_deflated_reproductive_frac",
    "rna_deflated_reproductive_and_thymus_frac",
    "rna_max_ntpm",
    "rna_thymus",
]
_PCT_FILTERS = [("80", 0.80), ("90", 0.90), ("95", 0.95), ("97", 0.97), ("98", 0.98), ("99", 0.99)]


def _recompute_protein_columns(seed: pd.DataFrame, normal_tissue: pd.DataFrame) -> None:
    """Recompute the four HPA IHC protein columns in-place from normal_tissue.

    Rule (reverse-engineered to 97.5% reproduction of the shipped columns; the
    residual are genuine HPA release differences):

    * a gene absent from the IHC table OR detecting protein in no tissue
      (Level not in Low/Medium/High) -> all four columns are ``"no data"``
      (the original pipeline collapsed "antibody present, nothing detected"
      into ``no data``);
    * otherwise ``protein_strict_expression`` is the sorted detected tissues,
      ``protein_reliability`` the gene's HPA reliability tier,
      ``protein_reproductive`` whether the detected tissues are all
      reproductive (+thymus), and ``protein_thymus`` whether thymus is detected.
    """
    nt = normal_tissue.copy()
    nt["gid"] = nt["Gene"].astype(str)
    nt["tl"] = nt["Tissue"].astype(str).str.strip().str.lower()
    detected = nt[nt["Level"].isin(_PROTEIN_DETECTED_LEVELS)]
    det_by_gene = detected.groupby("gid")["tl"].apply(lambda s: sorted(set(s)))
    rel_by_gene = nt.groupby("gid")["Reliability"].first()

    for idx, row in seed.iterrows():
        gid = str(row["Ensembl_Gene_ID"]).split(".")[0]
        det = None if gid in _CROSS_REACTIVE_IHC else det_by_gene.get(gid)
        if not det:  # absent/cross-reactive IHC, or detected nowhere -> no data
            seed.at[idx, "protein_strict_expression"] = "no data"
            seed.at[idx, "protein_reliability"] = "no data"
            seed.at[idx, "protein_reproductive"] = "no data"
            seed.at[idx, "protein_thymus"] = "no data"
            continue
        seed.at[idx, "protein_strict_expression"] = "; ".join(det)
        seed.at[idx, "protein_reliability"] = rel_by_gene[gid]
        # Capitalized to match the shipped table's bool-string representation.
        seed.at[idx, "protein_reproductive"] = str(set(det) <= _PROTEIN_REPRODUCTIVE)
        seed.at[idx, "protein_thymus"] = str("thymus" in det)


def _recompute_rna_columns(df: pd.DataFrame, consensus: pd.DataFrame, consensus_path: Path):
    """Return a copy of *df* with all RNA-derived columns recomputed from HPA."""
    core = CORE_REPRODUCTIVE_TISSUES
    core_thymus = frozenset(core | {"thymus"})
    ntpm_by_gene: dict[str, dict[str, float]] = {}
    for ensg, sub in consensus.groupby("Gene"):
        ntpm_by_gene[str(ensg)] = {
            t.strip().lower(): float(v) for t, v in zip(sub["Tissue"], sub["nTPM"])
        }

    seed = df.copy()
    missing_from_consensus = []
    for idx, row in seed.iterrows():
        ensg = str(row["Ensembl_Gene_ID"]).split(".")[0]
        ntpm = ntpm_by_gene.get(ensg)
        if not ntpm:
            missing_from_consensus.append(row["Symbol"])
            continue
        deflated = _fraction(ntpm, core, deflate=True)
        seed.at[idx, "rna_reproductive_frac"] = round(_fraction(ntpm, core, deflate=False), 4)
        seed.at[idx, "rna_reproductive_and_thymus_frac"] = round(
            _fraction(ntpm, core_thymus, deflate=False), 4
        )
        seed.at[idx, "rna_deflated_reproductive_frac"] = round(deflated, 4)
        seed.at[idx, "rna_deflated_reproductive_and_thymus_frac"] = round(
            _fraction(ntpm, core_thymus, deflate=True), 4
        )
        for pct, thr in _PCT_FILTERS:
            seed.at[idx, f"rna_{pct}_pct_filter"] = bool(deflated >= thr)
        seed.at[idx, "rna_max_ntpm"] = round(max(ntpm.values()), 1)
        seed.at[idx, "rna_thymus"] = bool(ntpm.get("thymus", 0.0) >= 1.0)

    # Per-tissue detail + safety max columns + axes, exactly as add_cta_gene does.
    seed = tiers.enrich_rna_per_tissue(seed, str(consensus_path))
    seed["rna_reproductive"] = seed["rna_somatic_detected_count"].fillna(0).astype(int).eq(0)
    seed = tiers.assign_all_axes(seed)  # preserves existing ms_restriction column

    seed["passes_filters"] = seed.apply(
        lambda r: passes_filters_rule(r, HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS["Missing"]), axis=1
    )
    seed["never_expressed"] = seed.apply(
        lambda r: never_expressed_rule(r, HPA_EXPRESSION_FLOOR_NTPM), axis=1
    )
    if "filtered" in seed.columns:  # backward-compat alias
        seed["filtered"] = seed["passes_filters"]
    return seed, missing_from_consensus


def _expressed_set(df: pd.DataFrame) -> set[str]:
    """CTA_gene_names()-equivalent: passes_filters and not never_expressed."""
    pf = df["passes_filters"].astype(bool)
    ne = df["never_expressed"].astype(str).str.lower() == "true"
    return set(df.loc[pf & ~ne, "Symbol"])


def _report(old: pd.DataFrame, new: pd.DataFrame, missing: list[str], columns: list[str]) -> None:
    def boolcol(df, c):
        return df.set_index("Symbol")[c].astype(str).str.lower().eq("true")

    print("=" * 78)
    print("REGENERATION DELTA REPORT  (current HPA consensus vs shipped table)")
    print("=" * 78)
    print(f"rows: {len(old)}   genes missing from current consensus: {missing or 'none'}")

    for col in ["passes_filters", "never_expressed"]:
        o, n = boolcol(old, col), boolcol(new, col)
        common = o.index.intersection(n.index)
        flips = common[o[common] != n[common]]
        gained = sorted(flips[(~o[flips]) & n[flips]])
        lost = sorted(flips[o[flips] & (~n[flips])])
        print(f"\n{col}: {len(flips)} flips   +{len(gained)} / -{len(lost)}")
        if gained:
            print(f"  now TRUE : {', '.join(gained)}")
        if lost:
            print(f"  now FALSE: {', '.join(lost)}")

    # Net effect on the expressed default set (what CTA_gene_names returns).
    eo, en = _expressed_set(old), _expressed_set(new)
    print(f"\nCTA_gene_names() expressed set: {len(eo)} -> {len(en)}")
    if en - eo:
        print(f"  ADDED  ({len(en - eo)}): {', '.join(sorted(en - eo))}")
    if eo - en:
        print(f"  REMOVED({len(eo - en)}): {', '.join(sorted(eo - en))}")

    # Per-tissue / numeric column churn.
    print("\nper-column value changes (cells differing):")
    o_idx, n_idx = old.set_index("Symbol"), new.set_index("Symbol")
    common = o_idx.index.intersection(n_idx.index)
    for col in columns:
        if col in ("Symbol",) or col not in n_idx.columns:
            continue
        oc = o_idx.loc[common, col].astype(str)
        nc = n_idx.loc[common, col].astype(str)
        diff = int((oc != nc).sum())
        if diff:
            print(f"  {col:42s} {diff:4d} changed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consensus", default=None, help="Local HPA consensus path (else download).")
    parser.add_argument(
        "--normal-tissue", default=None, help="Local HPA normal_tissue path (else download)."
    )
    parser.add_argument("--apply", action="store_true", help="Overwrite the bundled CSV in place.")
    args = parser.parse_args()

    consensus_path = rna_tissue_consensus_path(args.consensus)
    consensus = pd.read_csv(consensus_path, sep="\t")
    normal_tissue = pd.read_csv(normal_tissue_path(args.normal_tissue), sep="\t")
    old = pd.read_csv(CSV_PATH)
    columns = list(old.columns)

    working = old.copy()
    _recompute_protein_columns(working, normal_tissue)  # IHC first; RNA axes read it
    new, missing = _recompute_rna_columns(working, consensus, consensus_path)
    new = new[columns]

    # Write first, then diff against the round-tripped CSV so the report reflects
    # what actually lands on disk (avoids spurious in-memory dtype "changes",
    # e.g. bool True vs str "True", NaN vs "").
    dest = CSV_PATH if args.apply else CSV_PATH.with_suffix(".regen.csv")
    new.to_csv(dest, index=False)
    _report(old, pd.read_csv(dest), missing, columns)

    if args.apply:
        print(f"\n--apply: wrote {len(new)} rows to {CSV_PATH}")
    else:
        print(f"\n(dry run) wrote regenerated table to {dest}; shipped CSV untouched.")
        print("Re-run with --apply to overwrite the bundled table.")


if __name__ == "__main__":
    main()
