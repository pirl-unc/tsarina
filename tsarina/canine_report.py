"""Sourced canine evidence tables, figures and an offline HTML report."""

from __future__ import annotations

import html
import json
import re
from hashlib import sha256
from pathlib import Path

import pandas as pd

from .vaccine_report import _fasta, _json_value


def _figures(result, out):
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams["svg.hashsalt"] = "tsarina-canine"
    import matplotlib.pyplot as plt

    figures = []
    curves = result["source_tables"]["cumulative_coverage"]
    if not curves.empty:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
        for label, field, color in (
            ("Tumor RNA", "expression", "#287a5b"),
            ("RNA + DLA MS support", "joint", "#305ca6"),
        ):
            axes[0].plot(
                curves.length_aa, curves[f"{field}_1_lower"] * 100, "o-", label=label, color=color
            )
            axes[0].fill_between(
                curves.length_aa,
                curves[f"{field}_1_lower"] * 100,
                curves[f"{field}_1_upper"] * 100,
                alpha=0.15,
                color=color,
            )
        axes[0].set(
            xlabel="Final-order construct prefix (aa)",
            ylabel="Declared tumor/genotype frame (%)",
            ylim=(0, 100),
        )
        axes[0].legend(fontsize=8)
        axes[1].stackplot(
            curves.length_aa,
            curves.observed_monoallelic_peptides,
            curves.inferred_only_peptides,
            labels=["Observed monoallelic restriction", "Affinity-inferred restriction only"],
            colors=["#287a5b", "#b48335"],
            alpha=0.8,
        )
        axes[1].set(
            xlabel="Final-order construct prefix (aa)",
            ylabel="Distinct MS-observed peptide strings",
        )
        axes[1].legend(fontsize=8)
        fig.savefig(out / "coverage-and-MS.svg", metadata={"Date": None})
        plt.close(fig)
        figures.append("coverage-and-MS.svg")
    ligands = result["ligands"]
    if not ligands.empty:
        retained = ligands[ligands.assembled.eq(True)]
        counts = [retained[retained.allele.eq(a)].peptide.nunique() for a in result["alleles"]]
        fig, ax = plt.subplots(figsize=(8, max(3, 0.35 * len(counts))), constrained_layout=True)
        ax.barh(result["alleles"], counts, color="#305ca6")
        ax.set(
            xlabel="Distinct retained MS-observed peptides with DLA support",
            title="DLA evidence; allele bars are not population frequencies",
        )
        fig.savefig(out / "dla-evidence.svg", metadata={"Date": None})
        plt.close(fig)
        figures.append("dla-evidence.svg")
    funnel = result["funnel"]
    if not funnel.empty:
        names = list(funnel.name)
        fig, ax = plt.subplots(figsize=(10, max(3, 0.5 * len(names))), constrained_layout=True)
        for field, label, color in (
            ("raw_aa", "Full protein", "#ddd9cf"),
            ("normal_ms_filtered_aa", "After specificity and healthy MS", "#8299b9"),
            ("assembled_aa", "Retained in construct", "#287a5b"),
        ):
            ax.barh(names, funnel[field], label=label, color=color)
        ax.set(xlabel="Native residues (linkers and initiation excluded)")
        ax.legend(fontsize=8)
        fig.savefig(out / "sequence-funnel.svg", metadata={"Date": None})
        plt.close(fig)
        figures.append("sequence-funnel.svg")
        fig, axes = plt.subplots(
            len(names),
            1,
            figsize=(11, max(3, 2 * len(names))),
            squeeze=False,
            constrained_layout=True,
        )
        for ax, row in zip(axes[:, 0], funnel.itertuples(index=False)):
            protein = result["ranking"].set_index("proteoform_key").loc[row.proteoform_key]
            ax.broken_barh([(0, len(protein.sequence))], (0, 0.5), facecolors="#ddd9cf")
            if result["design"]:
                for layer in result["design"]["layers"]:
                    if layer.get("proteoform_key") == row.proteoform_key:
                        ax.broken_barh(
                            [(layer["native_start"], layer["native_end"] - layer["native_start"])],
                            (0, 0.5),
                            facecolors="#287a5b",
                        )
            maps = result["source_tables"]["protein_ms_map"]
            for hit in (
                maps[maps.proteoform_key.eq(row.proteoform_key)].to_dict("records")
                if not maps.empty
                else []
            ):
                healthy = hit["mapping"] == "verified_healthy_MS_8mer_overlap"
                cancer = hit["sample_health"] == "tumor"
                y = 1 if healthy else 2 if cancer else 3
                color = "#ab3f3b" if healthy else "#305ca6" if cancer else "#b48335"
                ax.broken_barh(
                    [(hit["start"], hit["end"] - hit["start"])], (y, 0.4), facecolors=color
                )
            ax.set(
                xlim=(0, len(protein.sequence)),
                yticks=[0.2, 1.2, 2.2, 3.2],
                yticklabels=["Retained / full", "Healthy MS 8-mer", "Tumor MS", "Other MS"],
                title=row.name,
                xlabel="Native residue (zero-based)",
            )
        fig.savefig(out / "native-MS-maps.svg", metadata={"Date": None})
        plt.close(fig)
        figures.append("native-MS-maps.svg")
    return figures


def canine_counts(result, bundle):
    design = result["design"]
    selected = result["ranking"][result["ranking"].selected]
    ligands = result["ligands"]
    retained = ligands[ligands.assembled.eq(True)] if not ligands.empty else ligands
    known = (
        set(retained.loc[retained.evidence_tier.eq("monoallelic_ms"), "peptide"])
        if not retained.empty
        else set()
    )
    all_peptides = set(retained.peptide) if not retained.empty else set()
    counts = {
        "source_occurrences": len(bundle["occurrences"]),
        "exact_protein_groups": len(result["ranking"]),
        "identified_cohort_dogs": int(result["ranking"].total_identified_dogs.iloc[0]),
        "cohort_samples": len(result["source_tables"]["cohort_sample_audit"]),
        "assembled_proteins": len(selected),
        "assembled_native_pieces": sum(layer["kind"] == "cta_segment" for layer in design["layers"])
        if design
        else 0,
        "construct_length_aa": design["length_aa"] if design else 0,
        "construct_length_nt": design["length_nt"] if design else 0,
        "MS_observed_peptides": len(all_peptides),
        "observed_monoallelic_peptides": len(known),
        "inferred_only_peptides": len(all_peptides - known),
        "MS_observations": retained.observation_id.nunique() if not retained.empty else 0,
        "peptide_DLA_pairs": len(retained[["peptide", "allele"]].drop_duplicates())
        if not retained.empty
        else 0,
        "panel_alleles": len(result["alleles"]),
    }
    counts.update(
        {
            "admitted_groups": int(result["ranking"].restriction_status.eq("admitted").sum()),
            "rejected_groups": int(result["ranking"].restriction_status.eq("rejected").sum()),
            "unknown_groups": int(result["ranking"].restriction_status.eq("unknown").sum()),
            "rejected_MS_observations": len(result["source_tables"]["ms_modality_rejected"]),
            "selection_enabled_alleles": int(
                result["source_tables"]["dla_capabilities"].selection_enabled.sum()
            )
            if not result["source_tables"]["dla_capabilities"].empty
            else 0,
        }
    )
    return counts


def write_canine_report(result, bundle, output_dir):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("Output directory already contains files")
    tables = {
        "ranking": result["ranking"],
        "ligands": result["ligands"],
        "funnel": result["funnel"],
        **result["source_tables"],
    }
    for name, frame in tables.items():
        frame.to_csv(out / f"{name}.csv", index=False)
    (out / "input-bundle.json").write_text(json.dumps(bundle, indent=2, allow_nan=False) + "\n")
    design = result["design"]
    if design:
        _fasta(out / "UNAPPROVED-antigen.fasta", "canine_exploratory_unassessed", design["protein"])
        _fasta(
            out / "UNAPPROVED-coding.fasta",
            result["config"]["vaccine_type"],
            design["coding_sequence"],
        )
        _fasta(
            out / "UNAPPROVED-construct.fasta",
            result["config"]["vaccine_type"],
            design["nucleotide_sequence"],
        )
        for name in ("layers", "junctions", "cleavage", "background_overlaps"):
            pd.DataFrame(design[name]).to_csv(out / f"{name}.csv", index=False)
    figures = _figures(result, out)
    counts = result["counts"]
    text = [
        "# Canine cancer-antigen evidence and exploratory design",
        "",
        f"Cohort: **{result['cohort']['name']}**. Status: **{result['status']}**.",
        "",
        f"RNA cohort: {result['cohort']['description']}. Population: {result['cohort']['population']}.",
        "",
        "Rank exact protein sequences by absolute RNA expression in independent untreated primary-tumor dogs, then retain contiguous regions with exact positive-MS observations and qualifying DLA affinity evidence. The construct is unassessed for canine presentation, processing, safety and efficacy.",
        "",
        f"Expression threshold: {result['expression_policy']['threshold']} TPM. Repeated biopsies count once: all measured biopsies must pass for the lower bound; any possibly positive or missing biopsy contributes to the upper bound. Missing RNA is unknown. These are evidence bounds, not confidence intervals or clinical coverage.",
        "",
        f"Allowed normal tissues: {', '.join(result['restriction_policy']['allowed_tissues'])}. Policy: {result['restriction_policy']['name']} ({result['restriction_policy']['version']}); definition: {result['restriction_policy']['definition']}. Somatic RNA counterevidence threshold: {result['restriction_policy']['somatic_tpm_threshold']} TPM. These labels refer to the supplied canine policy; the human HPA definitions are not used.",
        "",
        "All rejected and unknown translated source occurrences remain specificity background. Verified healthy-primary MS outside the allowed tissues supplies an 8-mer exclusion. Tumor/adjacent tissue, cell lines and unverified tissue labels do not establish healthy-normal counterevidence; absence of observations does not establish safety.",
        "",
        "MS evidence belongs to the exact observed peptide. Monoallelic restriction is separate from affinity-inferred multiallelic/untyped assignment. Human-host DLA transfection establishes a distinct experimental context, not endogenous canine presentation. Raw sample and host/source taxonomy facets remain in the tables.",
        "",
        f"Genotype frame: {result['genotype_population']['name']}. {result['genotype_population']['description']}. Joint coverage uses per-target expression and MS-qualified alleles on full supplied genotypes. Pairing modes: {', '.join(result['coverage']['pairing_modes']) or 'none'}. Unsupported genotype mass: {result['coverage']['unsupported_genotype_mass']:.3f}; missing mass: {result['coverage']['missing_mass']:.3f}. Unsupported genotype mass counts dogs with any unassessed allele; it can overlap demonstrated reach through another allele and is not an additional coverage category. These estimates do not imply worldwide or breed coverage. No HWE, human mortality prior or human p95 is used.",
        "",
        f"Affinity provider: {result['prediction_provider']}. Capability and model/sequence hashes: `dla_capabilities.csv`. Cleavage is unassessed unless an explicit callback was supplied; injected predictions retain that provenance and do not validate a model. Coding policy: `{result['config']['codon_species']}`; generic coding uses a deterministic standard genetic-code reverse translation without species optimization. UTR and polyA settings and total AA/nt limits are recorded in the manifest.",
        "",
    ]
    if result["limit"]:
        text.extend([f"Design limit: {result['limit']}.", ""])
    text.extend(
        [
            "## Counts",
            "",
            *[f"- {k.replace('_', ' ')}: {v}" for k, v in counts.items()],
            "",
            "## Data sources",
            "",
            *[
                f"- {key}: {s['url']} — version {s['version']}; asset SHA256 `{s['sha256']}`; license {s['license']}."
                for key, s in result["sources"].items()
            ],
            "",
            f"Input file SHA256: `{result['provenance']['bundle_sha256']}`. Reference identity: `{result['provenance']['reference_key']}`. The input copy has the same canonical JSON identity; its formatting may differ from the original file.",
            "",
            "## Tables and figures",
            "",
            *[f"- [{name}.csv]({name}.csv)" for name in tables],
            *[f"- [{f}]({f})" for f in figures],
            "",
            "See: https://github.com/pirl-unc/tsarina/issues/208; https://github.com/pirl-unc/oncoref/issues/571; https://github.com/pirl-unc/hitlist/issues/660; https://github.com/pirl-unc/hitlist/issues/661; https://github.com/openvax/mhctools/issues/544; https://github.com/openvax/mhcflurry/issues/490.",
        ]
    )
    (out / "report.md").write_text("\n".join(text) + "\n")

    def paragraph(t):
        value = html.escape(t)
        value = re.sub(r"\*\*(.*?)\*\*", r"<strong>\1</strong>", value)
        return "<p>" + re.sub(r"`([^`]+)`", r"<code>\1</code>", value) + "</p>"

    paragraphs = "".join(
        paragraph(t) for t in text[2:] if t and not t.startswith(("#", "-", "See:"))
    )
    display = (
        result["ranking"].drop(columns=["sequence"]).to_html(index=False, escape=True, border=0)
    )
    downloads = "".join(
        f'<li><a href="{html.escape(p.name, quote=True)}">{html.escape(p.name)}</a></li>'
        for p in sorted(out.iterdir())
        if p.suffix in {".csv", ".json", ".fasta", ".md"}
    )
    downloads += '<li><a href="manifest.json">manifest.json</a></li>'
    images = "".join(
        f'<figure><img src="{f}" alt="{f.removesuffix(".svg").replace("-", " ")}" /></figure>'
        for f in figures
    )
    (out / "index.html").write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Canine cancer-antigen evidence</title><style>body{font:16px/1.55 system-ui;max-width:1180px;margin:40px auto;padding:0 24px;color:#25312d}h1{line-height:1.2}table{border-collapse:collapse;font-size:13px}td,th{padding:7px;border-bottom:1px solid #ddd;text-align:left}.scroll{overflow:auto}img{width:100%;height:auto}figure{margin:32px 0}a{color:#245a91}code{overflow-wrap:anywhere}</style><h1>Canine cancer-antigen evidence and exploratory design</h1>'
        + paragraphs
        + "<h2>Counts</h2>"
        + pd.DataFrame(list(counts.items()), columns=["Measure", "Count"]).to_html(
            index=False, border=0
        )
        + images
        + '<h2>Exact protein groups</h2><div class="scroll">'
        + display
        + "</div><h2>Data sources</h2><ul>"
        + "".join(
            f'<li>{html.escape(k)}: <a href="{html.escape(s["url"], quote=True)}">{html.escape(s["url"])}</a>; version {html.escape(s["version"])}; SHA256 {s["sha256"]}</li>'
            for k, s in result["sources"].items()
        )
        + "</ul><h2>Downloads</h2><ul>"
        + downloads
        + "</ul></html>"
    )
    hashes = {p.name: sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    (out / "manifest.json").write_text(
        json.dumps(_json_value({**result, "artifact_sha256": hashes}), indent=2, allow_nan=False)
        + "\n"
    )
