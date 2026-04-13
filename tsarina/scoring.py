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

"""MHC-peptide binding prediction via topiary + mhctools.

Topiary wraps mhctools' binding predictors (MHCflurry, NetMHCpan, …) with a
uniform DataFrame interface.  This module exposes a thin tsarina-shaped API
over that: one column per kind of prediction (presentation_score,
presentation_percentile, affinity_nm), wide rather than long.

Typical usage::

    from tsarina.scoring import score_presentation
    from tsarina.alleles import get_panel

    scores = score_presentation(
        peptides=["SLYNTVATL", "GILGFVFTL"],
        alleles=get_panel("iedb27_ab"),
        predictor="mhcflurry",
    )

Requires the relevant mhctools backend for the chosen predictor.  MHCflurry:
``pip install mhcflurry && mhcflurry-downloads fetch``.  NetMHCpan requires
the external binary from DTU.

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

import warnings

import pandas as pd

AFFINITY_THRESHOLDS_NM: tuple[int, ...] = (50, 150, 500)
PRESENTATION_SCORE_THRESHOLDS: tuple[float, ...] = (0.90, 0.95)
PRESENTATION_PERCENTILE_THRESHOLDS: tuple[float, ...] = (0.25, 0.50, 1.0)


def _resolve_predictor_class(predictor: str):
    """Map a predictor key to its mhctools class."""
    key = predictor.lower().replace("-", "_")
    try:
        if key == "mhcflurry":
            from mhctools import MHCflurry

            return MHCflurry
        if key in {"netmhcpan", "netmhcpan_ba"}:
            from mhctools import NetMHCpan

            return NetMHCpan
        if key == "netmhcpan_el":
            from mhctools import NetMHCpanEL

            return NetMHCpanEL
    except ImportError as e:
        raise ImportError(
            f"Predictor '{predictor}' requires mhctools + the underlying tool. "
            f"Install mhctools (`pip install mhctools`) and the backend "
            f"(`pip install mhcflurry` for mhcflurry; NetMHCpan for netmhcpan*)."
        ) from e
    raise ValueError(
        f"Unknown predictor '{predictor}'. Supported: 'mhcflurry', 'netmhcpan', 'netmhcpan_el'."
    )


def _pivot_topiary(result: pd.DataFrame) -> pd.DataFrame:
    """Pivot topiary's long (peptide, allele, kind) output to one row per pair.

    Topiary emits one row per (peptide, allele, kind) with kind in
    {``pMHC_affinity``, ``pMHC_presentation``}.  tsarina's callers expect a
    single row per (peptide, allele) with separate columns for each kind.
    """
    if result.empty:
        return pd.DataFrame(
            columns=[
                "peptide",
                "allele",
                "presentation_score",
                "presentation_percentile",
                "affinity_nm",
            ]
        )

    presentation = result[result["kind"] == "pMHC_presentation"][
        ["peptide", "allele", "score", "percentile_rank"]
    ].rename(columns={"score": "presentation_score", "percentile_rank": "presentation_percentile"})
    affinity = result[result["kind"] == "pMHC_affinity"][["peptide", "allele", "value"]].rename(
        columns={"value": "affinity_nm"}
    )

    out = presentation.merge(affinity, on=["peptide", "allele"], how="outer")
    cols = ["peptide", "allele", "presentation_score", "presentation_percentile", "affinity_nm"]
    return out[[c for c in cols if c in out.columns]].reset_index(drop=True)


def score_presentation(
    peptides: list[str],
    alleles: list[str],
    predictor: str = "mhcflurry",
    peptide_lengths: list[int] | None = None,
) -> pd.DataFrame:
    """Score peptide-allele pairs for MHC presentation via topiary.

    Parameters
    ----------
    peptides
        Peptide sequences.  Each peptide is scored against every allele.
    alleles
        HLA allele names (e.g. ``["HLA-A*02:01"]``).  Use mhcgnomes-friendly
        formatting; topiary passes through to mhctools for normalization.
    predictor
        Which mhctools backend to use.  One of ``"mhcflurry"`` (default),
        ``"netmhcpan"``, ``"netmhcpan_el"``.
    peptide_lengths
        Lengths to request when the backend needs them (defaults to the
        distinct lengths present in ``peptides``).

    Returns
    -------
    pd.DataFrame
        Columns: ``peptide``, ``allele``, ``presentation_score``,
        ``presentation_percentile``, ``affinity_nm``.  A backend that does
        not emit a given kind leaves that column NaN.
    """
    try:
        from topiary import TopiaryPredictor
    except ImportError as e:
        raise ImportError(
            "Topiary is required for scoring. Install with: pip install topiary"
        ) from e

    predictor_cls = _resolve_predictor_class(predictor)

    if peptide_lengths is None:
        peptide_lengths = sorted({len(p) for p in peptides}) or [9]

    tp = TopiaryPredictor(
        models=[predictor_cls(alleles=alleles, default_peptide_lengths=peptide_lengths)],
        alleles=alleles,
    )
    name_to_peptide = {f"p{i}": pep for i, pep in enumerate(peptides)}
    result = tp.predict_from_named_peptides(name_to_peptide)
    return _pivot_topiary(result)


def score_affinity(
    peptides: list[str],
    alleles: list[str],
    predictor: str = "mhcflurry",
) -> pd.DataFrame:
    """Score peptide-allele pairs for MHC binding affinity only.

    Thin convenience wrapper around :func:`score_presentation` that
    returns just ``affinity_nm``.
    """
    df = score_presentation(peptides, alleles, predictor=predictor)
    return df[["peptide", "allele", "affinity_nm"]].reset_index(drop=True)


def score_netmhcpan(
    peptides: list[str],
    alleles: list[str],
    netmhcpan_path: str = "netMHCpan",
) -> pd.DataFrame:
    """Score peptide-allele pairs using NetMHCpan via subprocess.

    .. deprecated::
        Prefer :func:`score_presentation` with ``predictor="netmhcpan"``,
        which uses mhctools' NetMHCpan wrapper and handles normalization
        consistently with the MHCflurry path.  This subprocess shim will be
        removed in a future release.
    """
    warnings.warn(
        "score_netmhcpan is deprecated; use "
        "score_presentation(..., predictor='netmhcpan') instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    import shutil
    import subprocess
    import tempfile

    if not shutil.which(netmhcpan_path):
        raise FileNotFoundError(
            f"NetMHCpan not found at '{netmhcpan_path}'. "
            "Download from https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/"
        )

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
