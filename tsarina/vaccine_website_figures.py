"""Exportable scientific figures for the vaccine results website."""

from __future__ import annotations

from hashlib import sha256

import pandas as pd


def coverage_figures(data, result, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    from .vaccine_coverage import TIERS

    def save(fig, name):
        if result["provenance"].get("synthetic") or any(result["prediction_callbacks"].values()):
            fig.suptitle("Illustrative inputs / injected model predictions", fontsize=10)
        fig.tight_layout()
        fig.savefig(out / f"{name}.svg", bbox_inches="tight")
        fig.savefig(out / f"{name}.png", dpi=140, bbox_inches="tight")
        plt.close(fig)

    colors = ["#087c87", "#4469c8", "#8d78b5"]
    cumulative = pd.DataFrame(data["coverage"])
    rank = cumulative[cumulative.axis.eq("protein")]
    x = rank.step.to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    for key, label, color in zip(
        ("all", "typed", "measured"),
        ("All assignments", "Measured + typed inference", "Measured monoallelic only"),
        colors,
    ):
        axes[0].plot(x, rank[f"hla_{key}"] * 100, marker="o", label=label, color=color)
    axes[0].set(
        title="Estimated HLA carrier reach · CIWD known-allele proxy",
        ylabel="Individuals with ≥1 supported allele (%)",
        ylim=(0, 100),
    )
    axes[0].legend(fontsize=8)
    for key, label, color in zip(
        ("mortality", "incidence"), ("Mortality-weighted", "Incidence-weighted"), colors
    ):
        lower, upper = rank[f"{key}_lower"].to_numpy() * 100, rank[f"{key}_upper"].to_numpy() * 100
        axes[1].plot(x, lower, label=f"{label} lower bound", color=color)
        axes[1].fill_between(x, lower, upper, alpha=0.15, color=color)
        axes[1].plot(x, upper, linestyle=":", color=color)
    axes[1].set(
        title="Expression-union bounds · weighted by global burden",
        ylabel="Share of global reference burden (%)",
        ylim=(0, 100),
    )
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.set_xlabel("Proteins added in mortality-priority order")
        ax.set_xticks(x)
        ax.grid(alpha=0.15)
    save(fig, "cumulative-proteins")

    segment = cumulative[cumulative.axis.eq("segment")]
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].stackplot(
        segment.length_aa,
        *[segment[f"peptides_{tier}"] for tier in TIERS],
        labels=["Measured monoallelic", "Typed-sample inferred", "Untyped-panel inferred"],
        colors=colors,
        step="post",
    )
    axes[0].set(
        title="Distinct MS-observed peptides in final-order construct prefixes",
        ylabel="Distinct observed peptides",
    )
    axes[0].legend(fontsize=8, loc="upper left")
    axes[1].plot(segment.length_aa, segment.pmhc, marker=".", color=colors[1])
    axes[1].set(
        xlabel="Construct prefix length (aa; includes any preceding M/linkers)",
        ylabel="Distinct peptide-HLA pairs",
    )
    for ax in axes:
        ax.grid(alpha=0.15)
    save(fig, "ligands-vs-length")

    hla = pd.DataFrame(data["hla"])
    fig, axes = plt.subplots(2, 1, figsize=(15, 7), sharex=True)
    x = np.arange(len(hla))
    heights = hla.published_global_frequency.where(
        hla.global_source_label.str.startswith("CIWD", na=False)
    )
    axes[0].bar(x, heights * 100, color=colors[0])
    axes[0].set(
        ylabel="CIWD allele frequency (%)",
        title="Panel allele frequency and retained peptide evidence",
    )
    bottom = np.zeros(len(hla))
    for tier, color in zip(TIERS, colors):
        values = hla[tier].to_numpy()
        axes[1].bar(x, values, bottom=bottom, label=tier.replace("_", " "), color=color)
        bottom += values
    axes[1].set(ylabel="Distinct peptides per HLA allele")
    axes[1].set_xticks(
        x,
        [
            a.replace("HLA-", "") + (" †" if pd.isna(heights.iloc[i]) else "")
            for i, a in enumerate(hla.allele)
        ],
        rotation=90,
        fontsize=7,
    )
    axes[1].legend(fontsize=8)
    axes[1].set_xlabel("† Frequency unavailable in CIWD Table A2; omitted from carrier estimate")
    save(fig, "hla-frequency-evidence")

    intervals = result["specific_intervals"]
    for protein in data["proteins"]:
        fig, ax = plt.subplots(figsize=(11, 4.4))
        ax.barh(0, protein["raw_aa"], color="#e0e7ef", height=0.55)
        specific = intervals[intervals.proteoform_key.eq(protein["proteoform_key"])]
        for row in specific.itertuples(index=False):
            ax.barh(2, row.end - row.start, left=row.start, height=0.55, color="#91b8c7")
        before = result["source_tables"].get("cta_specific_before_normal_ms", intervals)
        for row in before[before.proteoform_key.eq(protein["proteoform_key"])].itertuples(
            index=False
        ):
            ax.barh(1, row.end - row.start, left=row.start, height=0.55, color="#b6c9d7")
        for row in protein["layers"]:
            ax.barh(
                3,
                row["native_end"] - row["native_start"],
                left=row["native_start"],
                height=0.55,
                color=colors[0],
            )
        tracks = {
            "cancer": (4, "#4469c8"),
            "healthy_nonreproductive": (5, "#bf584d"),
            "healthy_reproductive": (6, "#b58a32"),
            "tumor_adjacent": (7, "#bc783e"),
            "unknown_or_other": (8, "#8693a1"),
        }
        seen = set()
        for hit in protein["tissue_hits"]:
            k = (hit["source_group"], hit["start"], hit["end"])
            if k in seen:
                continue
            seen.add(k)
            y, color = tracks[hit["source_group"]]
            ax.barh(
                y,
                hit["end"] - hit["start"],
                left=hit["start"],
                height=0.55,
                color=color,
                alpha=0.55,
            )
        ax.set_yticks(
            range(9),
            [
                "Full native protein",
                "CTA-specific regions",
                "After normal-MS gate",
                "Chosen native segments",
                "Cancer MS",
                "Healthy nonreproductive MS",
                "Healthy reproductive MS",
                "Tumor-adjacent MS",
                "Unknown / other MS",
            ],
        )
        ax.set(
            xlim=(0, protein["raw_aa"]),
            xlabel="Native protein position (zero-based, half-open)",
            title=f"{protein['name']}: retained sequence and exact MS peptide source observations",
        )
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.15)
        save(fig, "protein-map-" + sha256(protein["proteoform_key"].encode()).hexdigest()[:12])
