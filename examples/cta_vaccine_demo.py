"""Reproducible SYNTHETIC example; no cancer/prediction values are evidence.

Run after installing tsarina[vaccine]:
    python examples/cta_vaccine_demo.py /tmp/cta-vaccine-demo
"""

import sys
from pathlib import Path

import pandas as pd

from tsarina.vaccine import design_vaccine
from tsarina.vaccine_construct import VaccineConfig
from tsarina.vaccine_inputs import Protein, VaccineInputs

ALLELES = ("HLA-A*02:01", "HLA-A*24:02")


def synthetic_affinity(peptides, alleles):
    return pd.DataFrame(
        [
            {"peptide": peptide, "allele": allele, "affinity_nm": 5000.0}
            for peptide in peptides
            for allele in alleles
        ]
    )


def synthetic_cleavage(sequences):
    return {sequence: [0.5] * len(sequence) for sequence in sequences}


def synthetic_presentation(peptides, alleles, predictor):
    return pd.DataFrame(
        [
            {
                "peptide": peptide,
                "allele": allele,
                "presentation_percentile": 0.2,
                "presentation_score": 0.9,
                "affinity_nm": 50.0,
            }
            for peptide in peptides
            for allele in alleles
        ]
    )


def main(output_dir):
    from unittest.mock import patch

    sequence = "MMMMMMMMACDEFGHIKLLLLLLLL"
    inputs = VaccineInputs(
        proteins=[
            Protein("demo1", "DEMO_A", "demo_p1", sequence),
            Protein("demo2", "DEMO_B", "demo_p2", sequence),
            Protein("demo3", "DEMO_C", "demo_p3", "MMMMMMMMTVWYACDEFGLLLLLLLL"),
            Protein("background", "NON_CTA", "background_p", "MMMMMMMMLLLLLLLL"),
        ],
        cta_gene_ids={"demo1", "demo2", "demo3"},
        gene_keys={"demo1": "DEMO_A/B", "demo2": "DEMO_A/B", "demo3": "demo3"},
        prevalence=pd.DataFrame(
            [
                {
                    "proteoform_key": key,
                    "cancer_code": code,
                    "prevalence_p95": prevalence,
                    "n_samples": 100,
                }
                for key, code, prevalence in [
                    ("DEMO_A/B", "DEMO_LUNG", 0.8),
                    ("DEMO_A/B", "DEMO_BREAST", 0.1),
                    ("demo3", "DEMO_LUNG", 0.2),
                    ("demo3", "DEMO_BREAST", 0.7),
                ]
            ]
        ),
        burden=pd.DataFrame(
            [
                {
                    "burden_category": "lung",
                    "world_mortality_pct": 20.0,
                    "world_incidence_pct": 12.0,
                    "source": "SYNTHETIC",
                },
                {
                    "burden_category": "breast",
                    "world_mortality_pct": 7.0,
                    "world_incidence_pct": 11.0,
                    "source": "SYNTHETIC",
                },
            ]
        ),
        cohorts={"lung": ["DEMO_LUNG"], "breast": ["DEMO_BREAST"]},
        ms_hits=pd.DataFrame(
            [
                {
                    "peptide": "ACDEFGHIK",
                    "mhc_restriction": ALLELES[0],
                    "is_monoallelic": True,
                    "pmid": "SYNTHETIC",
                    "assay_modality": "mass_spectrometry",
                },
                {
                    "peptide": "TVWYACDEF",
                    "mhc_restriction": ALLELES[1],
                    "is_monoallelic": True,
                    "pmid": "SYNTHETIC",
                    "assay_modality": "mass_spectrometry",
                },
            ]
        ),
        provenance={
            "synthetic": True,
            "warning": "All expression, burden, MS and predictions are illustrative fixtures",
        },
    )
    with patch("tsarina.scoring.score_presentation", synthetic_presentation):
        result = design_vaccine(
            VaccineConfig(top_k=2, alleles=ALLELES),
            inputs,
            output_dir=output_dir,
            affinity_fn=synthetic_affinity,
            cleavage_fn=synthetic_cleavage,
        )
    print(result["funnel"].to_string(index=False))
    print(f"SYNTHETIC audit: {Path(output_dir) / 'report.md'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/tmp/cta-vaccine-demo")
