"""Regenerate the gene-level MS safety overlay from Hitlist 2.x MS evidence.

``tsarina/data/gene-ms-safety-evidence.csv`` summarizes, per gene, where public
class I immunopeptidomics has observed peptides of that gene's reference
protein: cancer only, expected reproductive/thymic tissue, or healthy somatic
tissue. It has to be recomputed whenever the OncoRef target selection or the
indexed MS corpus changes.

The genes searched are every current OncoRef ``target_selection`` candidate
(``CTA`` plus ``Cancer Antigen only`` and ``Reproductive Antigen only``) plus
every gene the overlay already tracks, so a gene that leaves the selection
keeps evidence computed from the same corpus as everything else. Only genes
with at least one observed peptide get a row; the loader reports the rest as
``NO_MS_DATA``.

The sidecar ``gene-ms-safety-evidence.meta.json`` records which selection,
MS backend/corpus and sequence reference the rows describe. The test suite
compares it with the installed selection, because a member gene admitted after
the overlay was computed would otherwise read ``NO_MS_DATA`` without having
been searched.

Run against a current Hitlist 2.x observations index or verified EvidenceStore
package::

    python scripts/regenerate_gene_ms_safety_evidence.py
    python scripts/regenerate_gene_ms_safety_evidence.py --evidence-package /path/to/pkg
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

DEFAULT_OUTPUT_DIR = Path("tsarina/data")
OVERLAY_NAME = "gene-ms-safety-evidence.csv"
META_NAME = "gene-ms-safety-evidence.meta.json"
OVERLAY_COLUMNS = ("Ensembl_Gene_ID", "ms_restriction", "ms_healthy_somatic_tissues", "ms_pmids")


def compute_overlay(
    previous: pd.DataFrame | None = None,
    ensembl_release: int = 112,
    evidence_package: str | Path | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Return the overlay rows and the counts describing how they were derived."""
    from tsarina.gene_sets import CTA_gene_ids, CTA_target_candidate_gene_ids
    from tsarina.ms_evidence import load_public_ms_hits
    from tsarina.peptides import gene_peptides
    from tsarina.tiers import aggregate_gene_ms_safety

    member_ids = set(CTA_gene_ids())
    candidate_ids = set(CTA_target_candidate_gene_ids())
    tracked_ids = set() if previous is None else set(previous["Ensembl_Gene_ID"])
    searched_ids = sorted(member_ids | candidate_ids | tracked_ids)
    peptides = gene_peptides(searched_ids, ensembl_release=ensembl_release)
    hit_kwargs: dict = {
        "peptides": set(peptides["peptide"]),
        "mhc_class": "I",
        "classify_source": True,
        "drop_binding_assays": True,
    }
    if evidence_package is not None:
        hit_kwargs["evidence_package"] = evidence_package
    hits = load_public_ms_hits(**hit_kwargs)

    # Aggregate by gene ID rather than symbol: selection symbols and the
    # sequence release's gene names can differ, and some members have no symbol.
    gene_map = (
        peptides[["peptide", "gene_id"]].drop_duplicates().rename(columns={"gene_id": "gene_name"})
    )
    overlay = (
        aggregate_gene_ms_safety(hits, gene_map)
        .rename(columns={"gene_name": "Ensembl_Gene_ID"})
        .sort_values("Ensembl_Gene_ID")
        .reset_index(drop=True)
    )[list(OVERLAY_COLUMNS)]

    with_peptides = set(peptides["gene_id"])
    observed_members = overlay.loc[overlay["Ensembl_Gene_ID"].isin(member_ids)]
    observed_candidates = overlay.loc[overlay["Ensembl_Gene_ID"].isin(candidate_ids)]
    counts = {
        "genes_searched": len(searched_ids),
        "member_genes": len(member_ids),
        "member_genes_with_reference_peptides": len(member_ids & with_peptides),
        "member_genes_observed": len(observed_members),
        "candidate_genes": len(candidate_ids),
        "candidate_genes_with_reference_peptides": len(candidate_ids & with_peptides),
        "candidate_genes_observed": len(observed_candidates),
        "non_member_genes_observed": len(overlay) - len(observed_members),
        "reference_peptides_searched": int(peptides["peptide"].nunique()),
        "peptides_observed": int(hits["peptide"].nunique()) if not hits.empty else 0,
        "member_ms_restriction": (
            observed_members["ms_restriction"].value_counts().sort_index().to_dict()
        ),
    }
    return overlay, counts


def corpus_provenance(evidence_package: str | Path | None = None) -> dict:
    """Describe the Hitlist 2.x EvidenceStore or observations index used."""
    from tsarina.evidence_store import evidence_store_provenance, resolve_evidence_package

    resolved = resolve_evidence_package(evidence_package)
    if resolved is not None:
        return evidence_store_provenance(resolved)

    import hitlist
    import mhcgnomes
    from hitlist.observations import observations_path

    meta_path = observations_path().parent / "observations_meta.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    sources = {
        name: {key: source[key] for key in ("sha256", "size") if key in source}
        for name, source in meta.get("sources", {}).items()
        if name in ("iedb", "cedar")
    }
    return {
        "evidence_backend": "hitlist_observations",
        "hitlist_version": hitlist.__version__,
        "mhcgnomes_version": mhcgnomes.__version__,
        "observations_artifact_version": meta.get("artifact_version"),
        "observation_rows": meta.get("n_rows"),
        "observation_peptides": meta.get("n_peptides"),
        "raw_exports": sources,
    }


def build_meta(
    counts: dict,
    ensembl_release: int,
    evidence_package: str | Path | None = None,
) -> dict:
    import oncoref

    from tsarina.cta_selection import selection_provenance

    selection = selection_provenance(ensembl_release=ensembl_release)
    return {
        "oncoref_version": oncoref.__version__,
        "cta_selection_policy_ids": selection["cta_selection_policy_ids"],
        "cta_selection_reference_ids": selection["cta_selection_reference_ids"],
        "cta_selection_source_snapshot_ids": selection["cta_selection_source_snapshot_ids"],
        "cta_selection_release_ready": selection["cta_selection_release_ready"],
        "sequence_ensembl_release": ensembl_release,
        "mhc_class": "I",
        "mhc_species": "Homo sapiens",
        "corpus": corpus_provenance(evidence_package=evidence_package),
        "counts": counts,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--ensembl-release", type=int, default=112)
    parser.add_argument(
        "--evidence-package",
        type=Path,
        default=None,
        help="Optional verified Hitlist 2.x EvidenceStore package directory.",
    )
    args = parser.parse_args(argv)

    overlay_path = args.output_dir / OVERLAY_NAME
    previous = (
        pd.read_csv(overlay_path, dtype=str, keep_default_na=False)
        if overlay_path.exists()
        else None
    )
    overlay, counts = compute_overlay(
        previous,
        ensembl_release=args.ensembl_release,
        evidence_package=args.evidence_package,
    )
    meta = build_meta(counts, args.ensembl_release, evidence_package=args.evidence_package)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    overlay.to_csv(overlay_path, index=False)
    (args.output_dir / META_NAME).write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    print(json.dumps(meta["counts"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
