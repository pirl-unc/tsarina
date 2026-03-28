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

"""MHC-peptide binding prediction using MHCflurry.

Requires the ``mhcflurry`` package (not a perseo dependency -- install
separately with ``pip install mhcflurry`` and download models with
``mhcflurry-downloads fetch``).

Typical usage::

    from perseo.scoring import score_presentation
    from perseo.alleles import get_panel

    scores = score_presentation(
        peptides=["SLYNTVATL", "GILGFVFTL"],
        alleles=get_panel("iedb27_ab"),
    )

Thresholds commonly used for filtering:

.. list-table::
   :header-rows: 1

   * - Metric
     - Threshold
     - Meaning
   * - Affinity (nM)
     - <= 500
     - Weak binder
   * - Affinity (nM)
     - <= 150
     - Moderate binder
   * - Affinity (nM)
     - <= 50
     - Strong binder
   * - Presentation score
     - >= 0.90
     - Likely presented
   * - Presentation percentile
     - <= 1.0%
     - Top 1% of peptides for this allele
   * - Presentation percentile
     - <= 0.5%
     - Top 0.5% of peptides for this allele
"""

from __future__ import annotations

import pandas as pd

# Standard affinity and presentation thresholds used in CTA analysis
AFFINITY_THRESHOLDS_NM: tuple[int, ...] = (50, 150, 500)
PRESENTATION_SCORE_THRESHOLDS: tuple[float, ...] = (0.90, 0.95)
PRESENTATION_PERCENTILE_THRESHOLDS: tuple[float, ...] = (0.25, 0.50, 1.0)


def score_presentation(
    peptides: list[str],
    alleles: list[str],
    include_affinity: bool = True,
    include_processing: bool = False,
    n_flanks: list[str] | None = None,
    c_flanks: list[str] | None = None,
) -> pd.DataFrame:
    """Score peptide-allele pairs for MHC presentation using MHCflurry.

    Parameters
    ----------
    peptides
        List of peptide sequences.
    alleles
        List of HLA allele names (e.g. ``["HLA-A*02:01"]``).  Each
        peptide is scored against every allele.
    include_affinity
        If True (default), include ``affinity_nm`` and ``affinity_score``
        columns.
    include_processing
        If True, include ``processing_score`` column.  Requires flanking
        sequences for accurate processing prediction.
    n_flanks
        N-terminal flanking sequences (one per peptide).  Required if
        ``include_processing=True``.
    c_flanks
        C-terminal flanking sequences (one per peptide).  Required if
        ``include_processing=True``.

    Returns
    -------
    pd.DataFrame
        Columns: ``peptide``, ``allele``, ``presentation_score``,
        ``presentation_percentile``, and optionally ``affinity_nm``,
        ``affinity_score``, ``processing_score``.

    Raises
    ------
    ImportError
        If ``mhcflurry`` is not installed.
    """
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError as e:
        raise ImportError(
            "MHCflurry is required for scoring. Install with: pip install mhcflurry\n"
            "Then download models: mhcflurry-downloads fetch"
        ) from e

    predictor = Class1PresentationPredictor.load()

    # Build cross-product
    all_peptides = []
    all_alleles = []
    all_n_flanks = []
    all_c_flanks = []
    for i, pep in enumerate(peptides):
        for allele in alleles:
            all_peptides.append(pep)
            all_alleles.append(allele)
            if n_flanks is not None:
                all_n_flanks.append(n_flanks[i])
            if c_flanks is not None:
                all_c_flanks.append(c_flanks[i])

    predict_kwargs: dict = {
        "peptides": all_peptides,
        "alleles": all_alleles,
    }
    if n_flanks is not None:
        predict_kwargs["n_flanks"] = all_n_flanks
    if c_flanks is not None:
        predict_kwargs["c_flanks"] = all_c_flanks

    result_df = predictor.predict(**predict_kwargs)

    # Rename columns to standard names
    rename = {
        "presentation_score": "presentation_score",
        "presentation_percentile": "presentation_percentile",
    }
    if include_affinity:
        rename["affinity"] = "affinity_nm"

    cols_to_keep = ["peptide", "allele", "presentation_score", "presentation_percentile"]
    if include_affinity and "affinity" in result_df.columns:
        result_df = result_df.rename(columns={"affinity": "affinity_nm"})
        cols_to_keep.append("affinity_nm")
    if include_processing and "processing_score" in result_df.columns:
        cols_to_keep.append("processing_score")

    available = [c for c in cols_to_keep if c in result_df.columns]
    return result_df[available].reset_index(drop=True)


def score_affinity(
    peptides: list[str],
    alleles: list[str],
) -> pd.DataFrame:
    """Score peptide-allele pairs for MHC binding affinity only.

    Parameters
    ----------
    peptides
        List of peptide sequences.
    alleles
        List of HLA allele names.

    Returns
    -------
    pd.DataFrame
        Columns: ``peptide``, ``allele``, ``affinity_nm``,
        ``affinity_score`` (0-1, higher = stronger binding).

    Raises
    ------
    ImportError
        If ``mhcflurry`` is not installed.
    """
    try:
        from mhcflurry import Class1AffinityPredictor
    except ImportError as e:
        raise ImportError(
            "MHCflurry is required for scoring. Install with: pip install mhcflurry\n"
            "Then download models: mhcflurry-downloads fetch"
        ) from e

    predictor = Class1AffinityPredictor.load()

    all_peptides = []
    all_alleles = []
    for pep in peptides:
        for allele in alleles:
            all_peptides.append(pep)
            all_alleles.append(allele)

    result_df = predictor.predict_to_dataframe(
        peptides=all_peptides,
        alleles=all_alleles,
    )

    if "prediction" in result_df.columns:
        result_df = result_df.rename(columns={"prediction": "affinity_nm"})

    cols = [c for c in ["peptide", "allele", "affinity_nm"] if c in result_df.columns]
    return result_df[cols].reset_index(drop=True)


def score_netmhcpan(
    peptides: list[str],
    alleles: list[str],
    netmhcpan_path: str = "netMHCpan",
) -> pd.DataFrame:
    """Score peptide-allele pairs using NetMHCpan (external binary).

    Requires NetMHCpan 4.1+ installed and accessible at ``netmhcpan_path``.

    Parameters
    ----------
    peptides
        List of peptide sequences.
    alleles
        List of HLA allele names (e.g. ``["HLA-A02:01"]``).
        Note: NetMHCpan uses format without ``*`` (e.g. ``HLA-A02:01``).
    netmhcpan_path
        Path to the NetMHCpan binary (default ``"netMHCpan"``).

    Returns
    -------
    pd.DataFrame
        Columns: ``peptide``, ``allele``, ``affinity_nm``,
        ``rank_el`` (EL %Rank), ``rank_ba`` (BA %Rank).

    Raises
    ------
    FileNotFoundError
        If NetMHCpan binary is not found.
    RuntimeError
        If NetMHCpan fails.
    """
    import shutil
    import subprocess
    import tempfile

    if not shutil.which(netmhcpan_path):
        raise FileNotFoundError(
            f"NetMHCpan not found at '{netmhcpan_path}'. "
            "Download from https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/"
        )

    # NetMHCpan expects alleles without * (HLA-A02:01 not HLA-A*02:01)
    def _fmt_allele(a: str) -> str:
        return a.replace("*", "")

    rows: list[dict] = []
    for allele in alleles:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pep", delete=False) as f:
            for pep in peptides:
                f.write(pep + "\n")
            pep_file = f.name

        try:
            result = subprocess.run(
                [netmhcpan_path, "-p", pep_file, "-a", _fmt_allele(allele), "-BA"],
                capture_output=True,
                text=True,
                timeout=300,
            )
            if result.returncode != 0:
                raise RuntimeError(f"NetMHCpan failed: {result.stderr[:500]}")

            for line in result.stdout.splitlines():
                parts = line.strip().split()
                if len(parts) < 13 or not parts[0].isdigit():
                    continue
                try:
                    rows.append(
                        {
                            "peptide": parts[2],
                            "allele": allele,
                            "affinity_nm": float(parts[12]),
                            "rank_ba": float(parts[13]) if len(parts) > 13 else None,
                            "rank_el": float(parts[11]),
                        }
                    )
                except (ValueError, IndexError):
                    continue
        finally:
            import os

            os.unlink(pep_file)

    return pd.DataFrame(rows)
