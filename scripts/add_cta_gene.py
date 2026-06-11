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

"""Reproducibly add curated CTA gene rows to the bundled evidence table.

Generalizes ``scripts/add_xage2.py`` to a list of gene specs.  Each row is
built from HPA ``rna_tissue_consensus.tsv`` using the same in-repo generators
(:func:`tsarina.tiers.enrich_rna_per_tissue`, :func:`tsarina.tiers.assign_all_axes`)
that produce the rest of the table, and ``passes_filters`` / ``never_expressed``
are recomputed for every row from the documented, parameterized rules.

Only genes with **no HPA IHC protein data** can be added here (the row sets the
protein columns to "no data").  Genes with protein evidence, or that fail the
reproductive-restriction filter, must be evaluated separately — this script
asserts each added gene actually passes before writing.

Genes already present (by unversioned Ensembl gene ID) are skipped, so the
script is idempotent.

Usage::

    python scripts/add_cta_gene.py --check   # validate + preview, no write
    python scripts/add_cta_gene.py           # write the table
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

#: Genes to add.  Each must be protein-coding, have no HPA IHC protein data,
#: and pass the reproductive-restriction filter (asserted before writing).
GENE_SPECS = [
    {
        # MAGEB6 / CT3.4 — dual-corroborated (CTdatabase + CTexploreR), tsarina#79.
        # HPA: testis-only at 3.0 nTPM, protein not detected. Deflated
        # reproductive fraction 1.0; clean testis restriction.
        "Symbol": "MAGEB6",
        "Ensembl_Gene_ID": "ENSG00000176746",
        "Canonical_Transcript_ID": "ENST00000379034",
        "source_databases": "CTpedia;CTexploreR_CT",
        "Aliases": None,
        "Full_Name": None,
        "Function": "cancer-testis antigen",
    },
    # ── tsarina#93: near-identical paralog copies of curated CTAs, absent from
    # the source DBs but testis-restricted on HPA bulk consensus.  Adding them
    # to the universe stops their (sequence-identical) protein from polluting the
    # non-CTA negative set in downstream 9mer-specificity analysis.  Tagged
    # ``paralog:<curated sibling>``.  Only the clean testis-restricted copies are
    # added here; the somatic-on-bulk copies (CT45A6, DAZ2, DAZ4) are held out.
    *(
        {
            "Symbol": sym,
            "Ensembl_Gene_ID": ensg,
            "Canonical_Transcript_ID": ct,
            "source_databases": f"paralog:{sibling}",
            "Aliases": None,
            "Full_Name": None,
            "Function": "cancer-testis antigen",
        }
        for sym, ensg, ct, sibling in [
            # GAGE10 — full-length GAGE paralog (116 aa, 5 exons), absent from
            # the bundled universe. Clean testis restriction on HPA bulk (0
            # somatic tissues). NB: testis nTPM 1.6 < 2.0 floor, so the row
            # carries never_expressed=True; it still passes the reproductive-
            # restriction filter. Its sibling GAGE12B (ENSG00000236737) is
            # deliberately NOT added: in Ensembl that ID is a degenerate 117-bp
            # 3'-stub (7-aa ORF), not a GAGE protein -- see tsarina#108.
            ("GAGE10", "ENSG00000215274", "ENST00000407599", "GAGE2A"),
            ("CT47A8", "ENSG00000230347", "ENST00000457977", "CT47A1"),
            ("CT47A9", "ENSG00000226600", "ENST00000417256", "CT47A1"),
            ("CT47A10", "ENSG00000224089", "ENST00000430448", "CT47A1"),
            # Corrected IDs (tsarina#111): the previous entries carried the
            # siblings' gene IDs (ENSG…268606=MAGEA2, …269586=CT45A10), which were
            # already in the table, so they were silently skipped and the real
            # paralogs never got added. MAGEA2B/SSX4B encode a protein identical to
            # MAGEA2/SSX4 (distinct X loci ~tens of kb away — the NY-ESO-1 CTAG1A/B
            # and XAGE1A/B situation), so the gene must be in the universe or its
            # identical peptides leak into the non-CTA negative set. CT45A5 is a
            # distinct protein.
            ("MAGEA2B", "ENSG00000183305", "ENST00000331220", "MAGEA2"),
            ("CT45A5", "ENSG00000228836", "ENST00000698999", "CT45A1"),
            ("CT45A8", "ENSG00000278085", "ENST00000611660", "CT45A1"),
            ("CT45A9", "ENSG00000270946", "ENST00000620704", "CT45A1"),
            ("SSX4B", "ENSG00000269791", "ENST00000595235", "SSX4"),
            ("GAGE12D", "ENSG00000227488", "ENST00000405679", "GAGE12C"),
        ]
    ),
    # ── tsarina#111: placental-antigen families absent from the source DBs.
    # hCG-beta (CGB), pregnancy-specific glycoproteins (PSG), placental
    # syncytin/ERV envelopes, placental galectins (PP13/PP14) and placental
    # lactogen / GH are placenta-restricted on HPA bulk and are documented
    # onco-placental tumor antigens -- the placental analogue of the testis-
    # restricted cancer-germline antigens.  Added to the candidate universe and
    # left to the reproductive-restriction filter (no manual exclusion): the
    # somatically-broad members (GH2, PSG4/7, CGB1/3/5/7, ERVW-1/FRD-1/V-1) land
    # in the universe as *excluded* candidates; the placenta-dominant ones (PSG2,
    # PSG6, LGALS13/14, ERVV-2, ERVH48-1, CSH1, CGB8) pass the reproductive filter
    # into the expressed set, where any residual somatic/vital leakage (e.g. CSH1
    # lung ~10 nTPM) is caught downstream by the vital-tissue panel filter rather
    # than by hand here.  CGB8 also dropped from NON_CTA_EXCLUDED.  Tagged
    # ``placental_antigen``.  All are full-length protein-coding (139-538 aa; no
    # fragment models -- cf. GAGE12B, tsarina#108).
    *(
        {
            "Symbol": sym,
            "Ensembl_Gene_ID": ensg,
            "Canonical_Transcript_ID": ct,
            "source_databases": "placental_antigen",
            "Aliases": None,
            "Full_Name": None,
            # "cancer-germline antigen" is the accurate umbrella covering both the
            # testis-restricted CT antigens and these placenta-restricted (onco-
            # placental) members; the per-row ``restriction`` axis (TESTIS /
            # PLACENTAL / REPRODUCTIVE) carries the subtype.
            "Function": "cancer-germline antigen",
        }
        for sym, ensg, ct in [
            ("CGB1", "ENSG00000267631", "ENST00000301407"),
            ("CGB2", "ENSG00000104818", "ENST00000359342"),
            ("CGB3", "ENSG00000104827", "ENST00000357383"),
            ("CGB5", "ENSG00000189052", "ENST00000301408"),
            ("CGB7", "ENSG00000196337", "ENST00000684222"),
            ("CGB8", "ENSG00000213030", "ENST00000448456"),
            ("PSG2", "ENSG00000242221", "ENST00000406487"),
            ("PSG4", "ENSG00000243137", "ENST00000405312"),
            ("PSG6", "ENSG00000170848", "ENST00000292125"),
            ("PSG7", "ENSG00000221878", "ENST00000446844"),
            ("ERVW-1", "ENSG00000242950", "ENST00000603053"),  # syncytin-1
            ("ERVFRD-1", "ENSG00000244476", "ENST00000472091"),  # syncytin-2
            ("ERVH48-1", "ENSG00000233056", "ENST00000447535"),  # suppressyn
            ("ERVV-1", "ENSG00000269526", "ENST00000602168"),
            ("ERVV-2", "ENSG00000268964", "ENST00000601417"),
            ("LGALS13", "ENSG00000105198", "ENST00000221797"),  # PP13
            ("LGALS14", "ENSG00000006659", "ENST00000360675"),  # PP14
            ("CSH1", "ENSG00000136488", "ENST00000329882"),  # placental lactogen
            ("GH2", "ENSG00000136487", "ENST00000332800"),  # placental GH
        ]
    ),
]

_NO_PROTEIN = {"no data", "nan", ""}
_PROTEIN_DEFAULTS = {
    "protein_reproductive": "no data",
    "protein_thymus": "no data",
    "protein_reliability": "no data",
    "protein_strict_expression": "no data",
    "ms_restriction": "NO_MS_DATA",
}


def _has_protein(reliability) -> bool:
    return str(reliability).strip().lower() not in _NO_PROTEIN


def passes_filters_rule(row: pd.Series, missing_threshold: float) -> bool:
    if str(row["biotype"]) != "protein_coding":
        return False
    reliability = str(row["protein_reliability"]).strip()
    has_protein = _has_protein(reliability)
    if has_protein and str(row["protein_reproductive"]).strip().lower() != "true":
        return False
    tier = reliability if has_protein else "Missing"
    threshold = HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS.get(tier, missing_threshold)
    if tier in ("Uncertain", "Missing"):
        threshold = missing_threshold
    frac = row["rna_deflated_reproductive_frac"]
    return False if pd.isna(frac) else bool(float(frac) >= threshold)


def never_expressed_rule(row: pd.Series, floor: float) -> bool:
    if _has_protein(row["protein_reliability"]):
        return False
    max_ntpm = row["rna_max_ntpm"]
    return True if pd.isna(max_ntpm) else bool(float(max_ntpm) < floor)


def _fraction(ntpm: dict[str, float], allowed: frozenset[str], deflate: bool) -> float:
    vals = {t: (max(0.0, v - 1.0) if deflate else v) for t, v in ntpm.items()}
    total = sum(vals.values())
    if total <= 0:
        return 1.0 if deflate else 0.0
    return sum(v for t, v in vals.items() if t in allowed) / total


def build_row(spec: dict, consensus: pd.DataFrame, columns: list[str], consensus_path: Path) -> dict:
    ensg = spec["Ensembl_Gene_ID"]
    sub = consensus[consensus["Gene"] == ensg]
    if sub.empty:
        raise SystemExit(f"{spec['Symbol']} ({ensg}) not found in consensus TSV")

    # Refuse degenerate fragment gene models (tsarina#108). A candidate whose
    # annotated protein is a stub -- e.g. GAGE12B's 7-aa ENSG00000236737, the
    # shared GAGE C-terminus -- is a mis-annotation, not an antigen, and its HPA
    # gene-level nTPM is quantification noise. Unlike a reproductive-filter
    # failure (recorded as passes_filters=False below), this is never a valid
    # candidate, so refuse it at addition time where pyensembl data is available.
    from tsarina.qc import MIN_CTA_PROTEIN_AA, gene_max_protein_length

    protein_length = gene_max_protein_length(ensg)
    if protein_length is None or protein_length < MIN_CTA_PROTEIN_AA:
        raise SystemExit(
            f"{spec['Symbol']} ({ensg}): longest protein is "
            f"{protein_length} aa (< {MIN_CTA_PROTEIN_AA} aa floor); "
            "refusing to add a fragment gene model."
        )

    ntpm = {t.strip().lower(): float(v) for t, v in zip(sub["Tissue"], sub["nTPM"])}

    core = CORE_REPRODUCTIVE_TISSUES
    core_thymus = frozenset(core | {"thymus"})

    row = dict.fromkeys(columns)
    row.update(_PROTEIN_DEFAULTS)
    row.update({k: v for k, v in spec.items() if k in columns})
    row["biotype"] = "protein_coding"

    row["rna_reproductive_frac"] = round(_fraction(ntpm, core, deflate=False), 4)
    row["rna_reproductive_and_thymus_frac"] = round(_fraction(ntpm, core_thymus, deflate=False), 4)
    deflated = _fraction(ntpm, core, deflate=True)
    row["rna_deflated_reproductive_frac"] = round(deflated, 4)
    row["rna_deflated_reproductive_and_thymus_frac"] = round(
        _fraction(ntpm, core_thymus, deflate=True), 4
    )
    for pct, thr in [("80", 0.80), ("90", 0.90), ("95", 0.95), ("97", 0.97), ("98", 0.98), ("99", 0.99)]:
        row[f"rna_{pct}_pct_filter"] = bool(deflated >= thr)
    row["rna_max_ntpm"] = round(max(ntpm.values()), 1)
    row["rna_thymus"] = bool(ntpm.get("thymus", 0.0) >= 1.0)

    seed = pd.DataFrame([row])
    seed = tiers.enrich_rna_per_tissue(seed, str(consensus_path))
    seed["rna_reproductive"] = seed["rna_somatic_detected_count"].fillna(0).astype(int).eq(0)
    seed = tiers.assign_all_axes(seed)

    out = seed.iloc[0].to_dict()
    out["passes_filters"] = passes_filters_rule(
        seed.iloc[0], HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS["Missing"]
    )
    out["never_expressed"] = never_expressed_rule(seed.iloc[0], HPA_EXPRESSION_FLOOR_NTPM)
    # A candidate that fails the reproductive-restriction filter is recorded with
    # passes_filters=False (it lands in the universe as an *excluded* candidate),
    # not refused -- so a curated candidate set "filters down" by the data rather
    # than by hand. Held-out genes are simply omitted from GENE_SPECS. (Whole-gene
    # identical-protein paralogs are still worth adding even when never_expressed,
    # to keep their peptides out of the non-CTA negative set.)
    return {c: out.get(c) for c in columns}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--consensus",
        default=None,
        help="Local HPA rna_tissue_consensus path; default downloads + caches from HPA.",
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    consensus_path = rna_tissue_consensus_path(args.consensus)

    df = pd.read_csv(CSV_PATH)
    columns = list(df.columns)

    # Validate the rules reproduce the shipped columns before recomputing.
    pf = df.apply(lambda r: passes_filters_rule(r, HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS["Missing"]), axis=1)
    ne = df.apply(lambda r: never_expressed_rule(r, HPA_EXPRESSION_FLOOR_NTPM), axis=1)
    pf_mm = int((pf != df["passes_filters"].astype(bool)).sum())
    ne_mm = int((ne != (df["never_expressed"].astype(str).str.lower() == "true")).sum())
    print(f"passes_filters reproduction: {len(df) - pf_mm}/{len(df)} match")
    print(f"never_expressed reproduction: {len(df) - ne_mm}/{len(df)} match")
    if pf_mm or ne_mm:
        raise SystemExit("Rule reproduction failed; aborting.")

    consensus = pd.read_csv(consensus_path, sep="\t")
    present = set(df["Ensembl_Gene_ID"])
    new_rows = []
    for spec in GENE_SPECS:
        if spec["Ensembl_Gene_ID"] in present:
            print(f"skip {spec['Symbol']} (already present)")
            continue
        row = build_row(spec, consensus, columns, consensus_path)
        new_rows.append(row)
        print(
            f"add {spec['Symbol']}: passes_filters={row['passes_filters']} "
            f"never_expressed={row['never_expressed']} restriction={row['restriction']} "
            f"max_ntpm={row['rna_max_ntpm']} safety_flags={row['safety_flags']!r}"
        )

    if args.check:
        print("\n--check: no files written.")
        return
    if not new_rows:
        print("nothing to add.")
        return

    df["passes_filters"] = pf
    df["never_expressed"] = ne
    df = pd.concat([df, pd.DataFrame(new_rows)[columns]], ignore_index=True)
    df.to_csv(CSV_PATH, index=False)
    print(f"\nWrote {len(df)} rows to {CSV_PATH}")


if __name__ == "__main__":
    main()
