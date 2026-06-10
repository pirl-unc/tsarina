# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Bundled HPA cancer-expression prevalence features for CTA ranking."""

from __future__ import annotations

from functools import lru_cache
from os.path import dirname, join

import pandas as pd

_DATA_DIR = join(dirname(__file__), "data")
_RNA_PREVALENCE_CSV = join(_DATA_DIR, "hpa-cancer-rna-prevalence.csv")
_IHC_PREVALENCE_CSV = join(_DATA_DIR, "hpa-cancer-ihc-prevalence.csv")
_DEFAULT_RNA_THRESHOLD = 2.0
_DEFAULT_CANCER_TYPE_PREVALENCE_FLOOR = 0.05
_TCGA_COHORT = "TCGA"
_TCGA_DEFAULT_THRESHOLDS = (1.0, 5.0)


def _threshold_label(value: float) -> str:
    return f"{value:g}".replace(".", "_").replace("-", "minus_")


@lru_cache(maxsize=1)
def hpa_cancer_rna_prevalence() -> pd.DataFrame:
    """Return bundled HPA TCGA/validation cancer RNA prevalence by gene/cancer."""
    return pd.read_csv(_RNA_PREVALENCE_CSV)


@lru_cache(maxsize=1)
def hpa_cancer_ihc_prevalence() -> pd.DataFrame:
    """Return bundled HPA cancer IHC prevalence by gene/cancer."""
    return pd.read_csv(_IHC_PREVALENCE_CSV)


def cta_cancer_expression_features(
    *,
    rna_threshold: float = _DEFAULT_RNA_THRESHOLD,
    cancer_type_prevalence_floor: float = _DEFAULT_CANCER_TYPE_PREVALENCE_FLOOR,
    rna_prevalence: pd.DataFrame | None = None,
    ihc_prevalence: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Summarize bundled HPA cancer prevalence into one row per CTA gene.

    ``cancer_rna_prevalent_cancer_type_count`` counts normalized cancer types
    where at least ``cancer_type_prevalence_floor`` of samples have RNA above
    ``rna_threshold`` pTPM. ``tumor_prevalence_panel_score`` is a transparent
    ranking score that prioritizes cancer-type breadth, then sample-level RNA
    prevalence, with IHC detection only as a weak tie-breaker.
    """
    rna = hpa_cancer_rna_prevalence() if rna_prevalence is None else rna_prevalence.copy()
    ihc = hpa_cancer_ihc_prevalence() if ihc_prevalence is None else ihc_prevalence.copy()

    rna_features = _rna_features(
        rna,
        rna_threshold=rna_threshold,
        cancer_type_prevalence_floor=cancer_type_prevalence_floor,
    )
    ihc_features = _ihc_features(ihc, cancer_type_prevalence_floor=cancer_type_prevalence_floor)

    if rna_features.empty and ihc_features.empty:
        return rna_features  # empty frame with [gene_id, symbol]
    if rna_features.empty:
        out = ihc_features
    elif ihc_features.empty:
        out = rna_features
    else:
        out = rna_features.merge(ihc_features, on=["gene_id", "symbol"], how="left")

    defaults = {
        # RNA features are absent when there is no cancer RNA prevalence for
        # these genes; default them so tumor_prevalence_panel_score still
        # computes (the function's documented output) instead of KeyError-ing.
        "cancer_rna_prevalent_cancer_type_count": 0,
        "cancer_rna_sample_prevalence": 0.0,
        "cancer_rna_max_cancer_type_prevalence": 0.0,
        "cancer_ihc_total_samples": 0,
        "cancer_ihc_detected_samples": 0,
        "cancer_ihc_medium_high_samples": 0,
        "cancer_ihc_detected_sample_prevalence": 0.0,
        "cancer_ihc_medium_high_sample_prevalence": 0.0,
        "cancer_ihc_detected_cancer_type_count": 0,
        "cancer_ihc_medium_high_cancer_type_count": 0,
    }
    for column, default in defaults.items():
        if column not in out.columns:
            out[column] = default
        out[column] = out[column].fillna(default)

    out["tumor_prevalence_panel_score"] = (
        out["cancer_rna_prevalent_cancer_type_count"].astype(float)
        + out["cancer_rna_sample_prevalence"].astype(float)
        + 0.01 * out["cancer_rna_max_cancer_type_prevalence"].astype(float)
        + 0.001 * out["cancer_ihc_medium_high_cancer_type_count"].astype(float)
    )
    return out


def _rna_features(
    rna: pd.DataFrame,
    *,
    rna_threshold: float,
    cancer_type_prevalence_floor: float,
) -> pd.DataFrame:
    if rna.empty:
        return pd.DataFrame(columns=["gene_id", "symbol"])

    label = _threshold_label(rna_threshold)
    count_col = f"expressed_samples_ptpm_ge_{label}"
    prevalence_col = f"prevalence_ptpm_ge_{label}"
    required = {"gene_id", "symbol", "cancer_type", "samples", count_col}
    missing = required - set(rna.columns)
    if missing:
        raise ValueError(
            "Bundled HPA cancer RNA prevalence table does not contain "
            f"threshold {rna_threshold:g} pTPM column(s): {sorted(missing)}"
        )

    rna = rna.copy()
    rna[count_col] = pd.to_numeric(rna[count_col], errors="coerce").fillna(0)
    rna["samples"] = pd.to_numeric(rna["samples"], errors="coerce").fillna(0)
    by_type = (
        rna.groupby(["gene_id", "symbol", "cancer_type"], as_index=False)
        .agg(
            samples=("samples", "sum"),
            expressed_samples=(count_col, "sum"),
            mean_ptpm=("mean_ptpm", "mean"),
            max_ptpm=("max_ptpm", "max"),
        )
        .assign(
            cancer_type_prevalence=lambda df: (
                df["expressed_samples"] / df["samples"].where(df["samples"] != 0, pd.NA)
            )
        )
    )
    if prevalence_col in rna.columns:
        by_label = (
            rna.groupby(["gene_id", "symbol"], as_index=False)[prevalence_col]
            .max()
            .rename(columns={prevalence_col: "cancer_rna_max_cohort_prevalence"})
        )
    else:
        by_label = pd.DataFrame(columns=["gene_id", "symbol", "cancer_rna_max_cohort_prevalence"])

    features = (
        by_type.groupby(["gene_id", "symbol"], as_index=False)
        .agg(
            cancer_rna_total_samples=("samples", "sum"),
            cancer_rna_expressed_samples=("expressed_samples", "sum"),
            cancer_rna_cancer_type_count=("cancer_type", "nunique"),
            cancer_rna_prevalent_cancer_type_count=(
                "cancer_type_prevalence",
                lambda values: int((values >= cancer_type_prevalence_floor).sum()),
            ),
            cancer_rna_mean_cancer_type_prevalence=("cancer_type_prevalence", "mean"),
            cancer_rna_max_cancer_type_prevalence=("cancer_type_prevalence", "max"),
            cancer_rna_max_ptpm=("max_ptpm", "max"),
        )
        .assign(
            cancer_rna_sample_prevalence=lambda df: (
                df["cancer_rna_expressed_samples"]
                / df["cancer_rna_total_samples"].where(df["cancer_rna_total_samples"] != 0, pd.NA)
            ),
            cancer_rna_threshold_ptpm=float(rna_threshold),
            cancer_type_prevalence_floor=float(cancer_type_prevalence_floor),
        )
    )
    if not by_label.empty:
        features = features.merge(by_label, on=["gene_id", "symbol"], how="left")
    return features.fillna(
        {
            "cancer_rna_sample_prevalence": 0.0,
            "cancer_rna_mean_cancer_type_prevalence": 0.0,
            "cancer_rna_max_cancer_type_prevalence": 0.0,
            "cancer_rna_max_cohort_prevalence": 0.0,
        }
    )


def _ihc_features(
    ihc: pd.DataFrame,
    *,
    cancer_type_prevalence_floor: float,
) -> pd.DataFrame:
    if ihc.empty:
        return pd.DataFrame(columns=["gene_id", "symbol"])

    required = {
        "gene_id",
        "symbol",
        "cancer",
        "total",
        "detected",
        "medium_high",
        "prevalence_detected",
        "prevalence_medium_high",
    }
    missing = required - set(ihc.columns)
    if missing:
        raise ValueError(
            f"Bundled HPA cancer IHC prevalence table is missing column(s): {sorted(missing)}"
        )

    ihc = ihc.copy()
    for column in ["total", "detected", "medium_high"]:
        ihc[column] = pd.to_numeric(ihc[column], errors="coerce").fillna(0)
    by_gene = ihc.groupby(["gene_id", "symbol"], as_index=False).agg(
        cancer_ihc_total_samples=("total", "sum"),
        cancer_ihc_detected_samples=("detected", "sum"),
        cancer_ihc_medium_high_samples=("medium_high", "sum"),
        cancer_ihc_detected_cancer_type_count=(
            "prevalence_detected",
            lambda values: int((values >= cancer_type_prevalence_floor).sum()),
        ),
        cancer_ihc_medium_high_cancer_type_count=(
            "prevalence_medium_high",
            lambda values: int((values >= cancer_type_prevalence_floor).sum()),
        ),
    )
    by_gene["cancer_ihc_detected_sample_prevalence"] = by_gene[
        "cancer_ihc_detected_samples"
    ] / by_gene["cancer_ihc_total_samples"].where(by_gene["cancer_ihc_total_samples"] != 0, pd.NA)
    by_gene["cancer_ihc_medium_high_sample_prevalence"] = by_gene[
        "cancer_ihc_medium_high_samples"
    ] / by_gene["cancer_ihc_total_samples"].where(by_gene["cancer_ihc_total_samples"] != 0, pd.NA)
    return by_gene.fillna(
        {
            "cancer_ihc_detected_sample_prevalence": 0.0,
            "cancer_ihc_medium_high_sample_prevalence": 0.0,
        }
    )


def cta_tcga_expression_features(
    *,
    thresholds: tuple[float, ...] = _TCGA_DEFAULT_THRESHOLDS,
    rna_prevalence: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Per-gene TCGA pan-cancer RNA prevalence summary.

    Aggregates the TCGA-cohort rows of the bundled HPA cancer RNA prevalence
    table into one row per ``(gene_id, symbol)`` with prevalence at each
    requested pTPM threshold (defaults: ``1.0`` and ``5.0``), the originating
    TCGA cancer type that maximizes prevalence at the lowest threshold, and
    counts of distinct TCGA cancer types contributing samples.

    The underlying TCGA expression measurements are the cohort=``TCGA`` rows
    of HPA v23's ``cancer RNA-Seq`` tab, which republishes harmonized TCGA
    Pan-Cancer pTPM values. Sample counts and 21 distinct TCGA project codes
    (BRCA, LUAD, ...) are preserved.
    """
    rna = hpa_cancer_rna_prevalence() if rna_prevalence is None else rna_prevalence
    if rna.empty or "cohort" not in rna.columns:
        return pd.DataFrame(columns=["gene_id", "symbol"])

    tcga = rna[rna["cohort"] == _TCGA_COHORT].copy()
    if tcga.empty:
        return pd.DataFrame(columns=["gene_id", "symbol"])

    if not thresholds:
        raise ValueError("cta_tcga_expression_features requires at least one threshold")
    sorted_thresholds = tuple(sorted(thresholds))
    threshold_labels = {t: _threshold_label(t) for t in sorted_thresholds}

    required = {"gene_id", "symbol", "cancer_type", "samples"}
    missing = required - set(tcga.columns)
    if missing:
        raise ValueError(
            f"Bundled HPA cancer RNA prevalence table missing TCGA column(s): {sorted(missing)}"
        )
    tcga["samples"] = pd.to_numeric(tcga["samples"], errors="coerce").fillna(0)
    for threshold, label in threshold_labels.items():
        count_col = f"expressed_samples_ptpm_ge_{label}"
        if count_col not in tcga.columns:
            raise ValueError(
                f"Bundled HPA cancer RNA prevalence table missing column {count_col!r} "
                f"required for TCGA threshold {threshold:g} pTPM"
            )
        tcga[count_col] = pd.to_numeric(tcga[count_col], errors="coerce").fillna(0)

    per_type = tcga.groupby(["gene_id", "symbol", "cancer_type"], as_index=False, sort=False).agg(
        samples=("samples", "sum"),
        **{
            f"expressed_samples_ptpm_ge_{label}": (f"expressed_samples_ptpm_ge_{label}", "sum")
            for label in threshold_labels.values()
        },
        max_ptpm=("max_ptpm", "max"),
    )
    for label in threshold_labels.values():
        count_col = f"expressed_samples_ptpm_ge_{label}"
        prev_col = f"cancer_type_prevalence_tpm_ge_{label}"
        per_type[prev_col] = per_type[count_col] / per_type["samples"].where(
            per_type["samples"] != 0, pd.NA
        )

    primary_label = threshold_labels[sorted_thresholds[0]]
    primary_prev_col = f"cancer_type_prevalence_tpm_ge_{primary_label}"
    top_type_col = f"tcga_top_cancer_type_tpm_ge_{primary_label}"
    top_prev_col = f"tcga_top_cancer_type_prevalence_tpm_ge_{primary_label}"

    ranked = per_type.sort_values(
        ["gene_id", "symbol", primary_prev_col],
        ascending=[True, True, False],
        na_position="last",
    )
    tops = (
        ranked.drop_duplicates(["gene_id", "symbol"], keep="first")[
            ["gene_id", "symbol", "cancer_type", primary_prev_col]
        ]
        .rename(
            columns={
                "cancer_type": top_type_col,
                primary_prev_col: top_prev_col,
            }
        )
        .reset_index(drop=True)
    )
    tops[top_type_col] = tops[top_type_col].astype(str).fillna("")
    tops[top_prev_col] = pd.to_numeric(tops[top_prev_col], errors="coerce").fillna(0.0)

    base_aggs: dict[str, tuple[str, str]] = {
        "tcga_sample_count": ("samples", "sum"),
        "tcga_cancer_type_count": ("cancer_type", "nunique"),
        "tcga_max_ptpm": ("max_ptpm", "max"),
    }
    for label in threshold_labels.values():
        base_aggs[f"tcga_expressed_samples_tpm_ge_{label}"] = (
            f"expressed_samples_ptpm_ge_{label}",
            "sum",
        )

    summary = per_type.groupby(["gene_id", "symbol"], as_index=False, sort=False).agg(**base_aggs)
    for label in threshold_labels.values():
        prev_col = f"tcga_pan_prevalence_tpm_ge_{label}"
        summary[prev_col] = summary[f"tcga_expressed_samples_tpm_ge_{label}"] / summary[
            "tcga_sample_count"
        ].where(summary["tcga_sample_count"] != 0, pd.NA)
        summary[prev_col] = summary[prev_col].fillna(0.0)

    summary = summary.merge(tops, on=["gene_id", "symbol"], how="left")
    summary[top_type_col] = summary[top_type_col].fillna("")
    summary[top_prev_col] = summary[top_prev_col].fillna(0.0)
    return summary
