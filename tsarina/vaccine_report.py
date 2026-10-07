"""Tables, sourced method report, sequence FASTAs and scientific figures."""

from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import pandas as pd


def _json_value(value):
    if isinstance(value, pd.DataFrame):
        return json.loads(value.to_json(orient="records"))
    if isinstance(value, dict):
        return {k: _json_value(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(v) for v in value]
    return value


def _fasta(path, label, sequence):
    path.write_text(
        f">{label}\n" + "\n".join(sequence[i : i + 80] for i in range(0, len(sequence), 80)) + "\n"
    )


def _table(frame, columns, formatters=None):
    formatters = formatters or {}
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in frame[columns].itertuples(index=False, name=None):
        cells = []
        for col, value in zip(columns, row):
            if pd.isna(value):
                cells.append("missing")
            elif col in formatters:
                cells.append(formatters[col](value))
            else:
                cells.append(str(value).replace("|", "\\|"))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_vaccine_report(result, output_dir):
    """Write a reproducible evidence bundle; manifest hashes every output file."""
    from .regions import allele_frequency_audit

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    # Avoid mixing artifacts from an old construct with a failed new design.
    if (out / "manifest.json").exists():
        raise ValueError(f"Output directory already contains a vaccine design: {out}")
    for name in (
        "ranking",
        "cancer_cohorts",
        "cancer_summary",
        "specific_intervals",
        "ligands",
        "funnel",
    ):
        result[name].to_csv(out / f"{name}.csv", index=False)
    for name, frame in result["source_tables"].items():
        frame.to_csv(out / f"{name}.csv", index=False)
    allele_frequency_audit(result["alleles"]).to_csv(out / "hla_panel.csv", index=False)
    design = result["design"]
    if design:
        empty_columns = {
            "junctions": [
                "peptide",
                "start",
                "end",
                "boundaries",
                "inside_linker",
                "allele",
                "affinity_nm",
                "below_threshold",
            ],
            "cleavage": ["bond", "context", "right_kind", "cleavage_probability"],
            "excluded_segments": ["segment_id", "reason"],
        }
        for name in ("layers", "junctions", "cleavage", "excluded_segments", "search_history"):
            frame = pd.DataFrame(design[name])
            if frame.empty and name in empty_columns:
                frame = pd.DataFrame(columns=empty_columns[name])
            frame.to_csv(out / f"{name}.csv", index=False)
        _fasta(out / "protein.fasta", "cta_vaccine_protein", design["protein"])
        _fasta(
            out / "cds.fasta",
            f"cta_vaccine_cds_{result['config']['vaccine_type']}",
            design["cds_nt"],
        )
        _fasta(
            out / "full.fasta",
            f"cta_vaccine_full_{result['config']['vaccine_type']}",
            design["full_nt"],
        )
    selected = result["ranking"].query("selected")
    names = set(selected.proteoform_key)
    cancer = result["cancer_summary"][result["cancer_summary"].proteoform_key.isin(names)]
    lines = [
        "# Mortality-prioritized CTA vaccine design",
        "",
        f"CTA definition: **{result['config']['definition']}**. Status: **{result['status']}**.",
        "",
        "## Method and interpretation",
        "",
        f"Selection mode: **{result['config']['selection_mode']}**; "
        + (
            "no protein count cap (top_k is ignored). "
            if result["config"]["selection_mode"] == "budget"
            else f"requested {result['config']['top_k']} proteoforms. "
        )
        + f"Exclusion patterns: {result['config']['exclude_gene_patterns']}; exact exceptions: {result['config']['allow_genes']}. `selection_screen.csv` records excluded, unsupported, length-rejected and uninspected candidates. Exclusions do not alter CTA membership or the non-CTA background. Supported mode reserves one whole ligand-bearing piece per contributing target before allocating extra pieces. Budget mode greedily allocates whole native pieces by marginal mortality-weighted lower expression-union gain per aa, then incidence gain, new panel alleles and observed peptide evidence. It is a heuristic, not a patient-level coverage optimum.",
        "",
        f"MS evidence source: **{result['provenance']['ms_input_kind'].replace('_', ' ')}**. Observation hashes and any supplied source snapshots/provenance are in `manifest.json` and the copied source tables.",
        "",
        "MS modality is positively established by the structured assay method, curated MS-only supplementary source, or declared supplied MS modality when no method is reported. Explicit non-MS and unknown-method records are excluded; `rejected_ms_observations.csv` preserves the source rows and reasons. Nonbinding is not equivalent to MS (Hitlist #644).",
        "",
        f"MS support mode: **{result['config']['ms_support_mode']}**. Sample-affinity mode retains the exact observed peptide with every typed sample allele predicted below {result['config']['ms_affinity_nm']:g} nM; it does not require best-of-haplotype or presentation-percentile selection. Untyped sample support by panel prediction: {result['config']['allow_untyped_ms']}. Study-wide allele pools are not sample genotypes. `ms_assignments.csv`, when present, links each observation to its measured or predicted assignment; `hla_support_counts.csv` reports distinct retained peptide counts by allele and evidence tier. Presentation mode uses the tier-percentile policy below. Neither mode infers a nested unobserved peptide from a longer observed sequence.",
        "",
        "The p95 prevalence is the fraction of cohort samples in which the proteoform is in the top 5% of that sample's expression ranking. Identical-sequence gene TPMs are summed BEFORE ranking by OncoRef. This is not a 95th-percentile TPM threshold across patients.",
        "",
        "Each proteoform score sums `world mortality share * p95 prevalence` across distinct cancer categories. Prevalence is sample-count-weighted across the specified broad cohorts; each mortality share is counted once. This prioritization is additive across proteoforms and does not estimate distinct patients covered, clinical benefit, or preventable deaths. The observed cohort mixture is not worldwide patient prevalence. Missing measurements are visible, and incomplete scores must be compared cautiously. Global cancer incidence is reported separately; it is not multiplied into the ranking score.",
        "",
        "Strict RNA numerator = OncoRef core reproductive tissues (testis, ovary, placenta). Loose adds cervix, endometrium, epididymis, fallopian tube, prostate, seminal vesicle and vagina. Neither RNA numerator adds breast or thymus. Thymus remains in the default fraction denominator but is excluded from somatic maxima and the reproductive protein flag; protein annotations use broader reproductive tissue conventions. These are curated RNA/protein evidence gates with exceptions, not absolute absence rules. Each definition has its own non-CTA background. Normal-tissue restriction does not establish target safety.",
        "",
        f"Normal-MS policy: **{result['config'].get('normal_ms_policy', 'audit')}**. When exclusion is enabled, primary donor-resolved nonmalignant Atlas HLA-I observations in heart, brain and lung blacklist peptides at {result['config'].get('normal_ms_min_donors', 1)} distinct donor(s). All residues covered by their 8-mers are removed after non-CTA subtraction. `normal_ms_summary.csv` and `normal_ms_audit.csv` preserve donor counts and qualification/exclusion reasons. These tissues are outside both CTA scopes. Other tissue observations remain warnings. Atlas autopsy donors had no diagnosed malignancy but could have other disease; nonmalignant is not disease-free. Dataset completeness and lack of observed presentation do not establish safety.",
        "",
        "All translated non-CTA coding isoforms are screened. Same-symbol HSCHR alternate-haplotype annotations inherit their curated primary gene's CTA identity, with exact resolutions in `background_cta_aliases.csv`; this does not change expression keys or admit independent non-CTA loci. Every residue covered by an independent non-CTA 8-mer is removed, including a full identical protein from another non-CTA gene. Native coordinates are zero-based, half-open. Each remaining contiguous interval must contain a qualifying panel MS ligand. In presentation mode, monoallelic evidence uses presentation percentile ≤2; sample-allele inference ≤1; unrestricted MS plus predicted assignment ≤0.5. Allele inference is labeled separately from measured restriction. Every qualifying ligand and repeated native occurrence is retained in the evidence tables.",
        "",
        "Padding trims terminal context only. Construct search compares complete sequences under the configured beam/round budget, including order, independent N/C padding and linker changes. Its lexicographic objective minimizes the number of peptide/allele junction predictions below the affinity cutoff, then their log binding burden, then maximizes mean predicted boundary cleavage, then minimizes linker length and retains context. The heuristic does not guarantee a global optimum. Final audit includes the initiating methionine, linker-internal windows and windows crossing multiple boundaries.",
        "",
        "Pepsickle predicts proteasomal cleavage with the human-only in-vivo model and eight residues of context on each side when available; a probability is a model output, not proof of cleavage. MHC binding does not establish presentation or immunogenicity. HLA frequencies are regional proxy evidence from the existing Tsarina panel audit, not a guarantee of population coverage.",
        "",
        "The coding sequence uses Vaxrank's DnaChisel codon optimization and is checked for exact translation, including initiation and stop. Total nucleotide length includes UTRs, CDS/stop and polyA. DNA outputs describe an insert and do not add a plasmid backbone or promoter. RNA FASTA uses U and does not encode a chemical cap or modified nucleosides.",
        "",
        "## Selected proteoforms",
        "",
        _table(
            selected,
            [
                "rank",
                "name",
                "gene_ids",
                "length_aa",
                "mortality_weighted_score",
                "complete_measurement",
            ],
            {"mortality_weighted_score": lambda v: f"{v:.6f}"},
        ),
        "",
        "## Cancer contributions",
        "",
        "Shares are percentages; p95 prevalence is a fraction. Absolute counts remain missing where OncoRef has no source counts. `cancer_cohorts.csv` includes each measured cohort and denominator; `cancer_summary.csv` includes every candidate and cancer category.",
        "",
        _table(
            cancer,
            [
                "name",
                "burden_category",
                "cancer_codes",
                "world_incidence_pct",
                "world_mortality_pct",
                "prevalence_p95",
                "n_samples",
                "complete_measurement",
                "score_contribution",
            ],
            {"prevalence_p95": lambda v: f"{v:.4f}", "score_contribution": lambda v: f"{v:.6f}"},
        ),
        "",
        "## Per-protein retention funnel",
        "",
        "MS-supported length includes the full specific interval before end trimming. Padded length uses maximum configured padding; assembled length reflects optimization and limits. Piece counts count native intervals, not unique peptide sequences.",
        "",
        _table(
            result["funnel"],
            [
                "name",
                "status",
                "raw_aa",
                "specific_aa",
                "specific_pieces",
                "normal_ms_filtered_aa",
                "normal_ms_filtered_pieces",
                "ms_supported_aa",
                "ms_supported_pieces",
                "padded_aa",
                "assembled_aa",
                "assembled_pieces",
                "assembled_ms_ligand_count",
                "assembled_pmhc_count",
                "retained_fraction",
                "assembled_panel_alleles",
            ],
            {"retained_fraction": lambda v: f"{v:.1%}"},
        ),
        "",
    ]
    if result["provenance"].get("synthetic"):
        lines[2:2] = [
            "**SYNTHETIC EXAMPLE: these cancer, expression, MS and prediction values are illustrative and must not be used as evidence.**",
            "",
        ]
    if design:
        bad = [r for r in design["junctions"] if r["below_threshold"]]
        lines.extend(
            [
                "## Construct and junction audit",
                "",
                f"Protein: **{design['length_aa']} aa**. Full {result['config']['vaccine_type'].upper()}: **{design['length_nt']} nt**. Junction affinity cutoff: **{result['config']['junction_affinity_nm']:g} nM**. Unresolved peptide/allele/window predictions: **{len(bad)}**.",
                "",
                "All predictions, including nonbinding windows, are in `junctions.csv`. Boundary cleavage probabilities are in `cleavage.csv`. `layers.csv` maps native source coordinates to assembled coordinates; `excluded_segments.csv` records whole-piece exclusions. A design with unresolved binders requires review; `--require-clean-junctions` makes that condition fail the command after writing the audit.",
                "",
            ]
        )
    lines.extend(
        [
            "## Figures and provenance",
            "",
            "![Cancer prevalence and mortality](cancer-priorities.svg)",
            "",
            "![Sequence retention](sequence-funnel.svg)",
            "",
        ]
    )
    if design:
        lines.extend(["![Construct layers](construct-map.svg)", ""])
    lines.extend(
        [
            "Exact configuration, dependency versions, source table hashes, sequence hashes, search history and predictor identity are in `manifest.json`. The copied reference tables retain upstream source anchors and provenance notes. Synthetic input/prediction callbacks are explicitly labeled.",
            "",
            "The friendly interactive report is in `website/index.html`. Serve this directory over HTTP to view it. Its downloads include cumulative cancer-expression union bounds, locus-aware HLA carrier proxies, MS counts versus final-order construct length and full native-protein tissue MS maps. Tissue observations are an audit layer; retained normal-tissue overlaps require review. Missing observations do not establish safety. Saved manifests can be rendered with `tsarina vaccine-report` without rerunning models.",
            "",
            "Scientific references: [OncoRef mortality report](https://github.com/pirl-unc/oncoref/blob/main/scripts/cta_mortality_coverage_report.py); [Vaxrank construction](https://github.com/openvax/vaxrank); [MHCflurry 2.0](https://pubmed.ncbi.nlm.nih.gov/32711842/); [Pepsickle](https://pubmed.ncbi.nlm.nih.gov/34478497/). OncoRef reference limitations: [mortality source/count refresh #542](https://github.com/pirl-unc/oncoref/issues/542), [category scope #543](https://github.com/pirl-unc/oncoref/issues/543).",
            "",
        ]
    )
    (out / "report.md").write_text("\n".join(lines))
    _figures(result, out)
    if design:
        from .vaccine_website import write_vaccine_website

        write_vaccine_website(
            {result["config"]["definition"]: result},
            out / "website",
            reports={result["config"]["definition"]: out},
        )
    hashes = {
        p.relative_to(out).as_posix(): sha256(p.read_bytes()).hexdigest()
        for p in sorted(out.rglob("*"))
        if p.is_file()
    }
    manifest = _json_value({**result, "artifact_sha256": hashes})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")


def _figures(result, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    selected = result["ranking"].query("selected")
    if selected.empty:
        return
    data = result["cancer_summary"]
    categories = list(
        dict.fromkeys(data.sort_values("world_mortality_pct", ascending=False).burden_category)
    )
    keys = selected.proteoform_key.tolist()
    matrix = data.pivot(
        index="proteoform_key", columns="burden_category", values="prevalence_p95"
    ).reindex(index=keys, columns=categories)
    fig, ax = plt.subplots(figsize=(max(10, len(categories) * 0.65), max(4, len(keys) * 0.45)))
    image = ax.imshow(matrix.values * 100, vmin=0, vmax=100, cmap="YlGnBu", aspect="auto")
    ax.set_yticks(range(len(keys)), selected.name)
    ref = data.drop_duplicates("burden_category").set_index("burden_category")
    ax.set_xticks(
        range(len(categories)),
        [f"{c}\n{ref.loc[c, 'world_mortality_pct']:g}% deaths" for c in categories],
        rotation=65,
        ha="right",
    )
    for y in range(len(keys)):
        for x in range(len(categories)):
            value = matrix.values[y, x]
            label = "NA" if np.isnan(value) else f"{100 * value:.0f}"
            ax.text(
                x,
                y,
                label,
                ha="center",
                va="center",
                fontsize=8,
                color="white" if value > 0.55 else "black",
            )
    ax.set_title("Observed p95 patient prevalence; mortality shares label cancer categories")
    fig.colorbar(image, ax=ax, label="Cohort samples in top 5% expression (%)")
    fig.tight_layout()
    if result["provenance"].get("synthetic"):
        fig.suptitle("SYNTHETIC EXAMPLE - illustrative values", y=1.08)
    fig.savefig(out / "cancer-priorities.svg", bbox_inches="tight")
    fig.savefig(out / "cancer-priorities.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    funnel = result["funnel"]
    stages = ["raw_aa", "specific_aa", "ms_supported_aa", "padded_aa", "assembled_aa"]
    if "normal_ms_filtered_aa" in funnel:
        stages.insert(2, "normal_ms_filtered_aa")
    fig, ax = plt.subplots(figsize=(10, max(3, len(funnel) * 0.65)))
    yy = np.arange(len(funnel))
    for i, stage in enumerate(stages):
        ax.barh(
            yy + (i - (len(stages) - 1) / 2) * 0.13,
            funnel[stage],
            height=0.12,
            label=stage.replace("_aa", "").replace("_", " "),
        )
    ax.set_yticks(yy, funnel.name)
    ax.invert_yaxis()
    ax.set_xlabel("Retained native amino acids (piece counts in funnel.csv)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    ax.set_title("Per-proteoform sequence retention")
    fig.tight_layout()
    if result["provenance"].get("synthetic"):
        fig.suptitle("SYNTHETIC EXAMPLE - illustrative values", y=1.08)
    fig.savefig(out / "sequence-funnel.svg", bbox_inches="tight")
    fig.savefig(out / "sequence-funnel.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    if result["design"]:
        design = result["design"]
        native = [layer for layer in design["layers"] if layer["kind"] == "cta_segment"]
        fig, ax = plt.subplots(figsize=(12, max(3, 1.5 + 0.4 * len(native))))
        colors = plt.get_cmap("tab20")
        keys = list(dict.fromkeys(layer["proteoform_key"] for layer in native))
        palette = {key: colors(i % 20) for i, key in enumerate(keys)}
        for layer in design["layers"]:
            ax.barh(
                0,
                layer["end_aa"] - layer["start_aa"],
                left=layer["start_aa"],
                color=palette.get(layer.get("proteoform_key"), "#777777"),
            )
        labels = ["Full antigen (added M/linkers in gray)"]
        for i, layer in enumerate(native, 1):
            ax.barh(
                i,
                layer["end_aa"] - layer["start_aa"],
                left=layer["start_aa"],
                color=palette[layer["proteoform_key"]],
            )
            labels.append(f"{i}. {layer['name']}  [{layer['native_start']}:{layer['native_end']}]")
        ax.set_xlim(0, design["length_aa"])
        ax.set_yticks(range(len(labels)), labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Construct position (aa); labels show native half-open intervals")
        ax.set_title(f"Single antigen: {design['length_aa']} aa / {design['length_nt']} total nt")
        fig.tight_layout()
        if result["provenance"].get("synthetic"):
            fig.suptitle("SYNTHETIC EXAMPLE - illustrative values", y=1.08)
        fig.savefig(out / "construct-map.svg", bbox_inches="tight")
        fig.savefig(out / "construct-map.png", dpi=160, bbox_inches="tight")
        plt.close(fig)
