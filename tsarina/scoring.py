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

"""MHC-peptide binding prediction.

MHCflurry is handled directly so tsarina can request only the values it uses
and batch peptide-allele presentation predictions efficiently. Other backends
still use topiary + mhctools' uniform DataFrame interface. This module exposes
a thin tsarina-shaped API over both paths: one column per kind of prediction
(``presentation_score``, ``presentation_percentile``, ``affinity_nm``), wide
rather than long.

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

.. note::
   A legacy ``score_netmhcpan`` subprocess shim lived here through v0.5.x.
   It was deleted in v0.6.0 after a side-by-side validation (HLA-A\\*02:01,
   5 known binders + 3 non-binder controls) showed its hand-rolled stdout
   parser returned the wrong output column on NetMHCpan 4.2c — affinities
   were off by ~1000x and qualitatively inverted (strong binders looked
   sub-nanomolar, polyA ranked like a moderate binder).  Topiary/mhctools'
   NetMHCpan adapter correctly detects the 4.2c output format, so
   ``score_presentation(..., predictor="netmhcpan")`` is the only path now.

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

import sys

import pandas as pd

AFFINITY_THRESHOLDS_NM: tuple[int, ...] = (50, 150, 500)
PRESENTATION_SCORE_THRESHOLDS: tuple[float, ...] = (0.90, 0.95)
PRESENTATION_PERCENTILE_THRESHOLDS: tuple[float, ...] = (0.25, 0.50, 1.0)

_SCORE_COLUMNS = [
    "peptide",
    "allele",
    "presentation_score",
    "presentation_percentile",
    "affinity_nm",
]
_MHCFLURRY_PRESENTATION_PREDICTOR: object | None = None


def _predictor_key(predictor: str) -> str:
    return predictor.lower().replace("-", "_")


def _empty_scores() -> pd.DataFrame:
    return pd.DataFrame(columns=_SCORE_COLUMNS)


def _unique_in_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def _load_mhcflurry_presentation_predictor():
    global _MHCFLURRY_PRESENTATION_PREDICTOR

    if _MHCFLURRY_PRESENTATION_PREDICTOR is None:
        try:
            from mhcflurry import Class1PresentationPredictor
        except ImportError as e:
            raise ImportError(
                "MHCflurry is required for predictor='mhcflurry' but is not "
                f"importable by the Python running tsarina ({sys.executable}). "
                "Install mhcflurry into that interpreter, or run tsarina with "
                "the environment where mhcflurry is installed."
            ) from e
        _MHCFLURRY_PRESENTATION_PREDICTOR = Class1PresentationPredictor.load()
    return _MHCFLURRY_PRESENTATION_PREDICTOR


def _score_mhcflurry(peptides: list[str], alleles: list[str]) -> pd.DataFrame:
    """Score with MHCflurry without requesting unused affinity percentiles.

    MHCflurry presentation percentiles are required by tsarina's ranking and
    tiering logic. Affinity percentile ranks are not used, and some supported
    alleles lack that separate calibration curve, so requesting them can fail
    even when raw affinity and presentation predictions are available.
    """
    unique_peptides = _unique_in_order(peptides)
    unique_alleles = _unique_in_order(alleles)
    if not unique_peptides or not unique_alleles:
        return _empty_scores()

    predictor = _load_mhcflurry_presentation_predictor()
    allele_samples = {allele: [allele] for allele in unique_alleles}
    result = predictor.predict(
        peptides=unique_peptides,
        alleles=allele_samples,
        include_affinity_percentile=False,
        verbose=0,
    )
    if result.empty:
        return _empty_scores()

    out = pd.DataFrame(
        {
            "peptide": result["peptide"],
            "allele": result["best_allele"] if "best_allele" in result else result["sample_name"],
            "presentation_score": result.get("presentation_score", pd.NA),
            "presentation_percentile": result.get("presentation_percentile", pd.NA),
            "affinity_nm": result.get("affinity", pd.NA),
        }
    )
    return out[_SCORE_COLUMNS].reset_index(drop=True)


def _resolve_predictor_class(predictor: str):
    """Map a predictor key to its mhctools class."""
    key = _predictor_key(predictor)
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
        return _empty_scores()

    presentation = result[result["kind"] == "pMHC_presentation"][
        ["peptide", "allele", "score", "percentile_rank"]
    ].rename(columns={"score": "presentation_score", "percentile_rank": "presentation_percentile"})
    affinity = result[result["kind"] == "pMHC_affinity"][["peptide", "allele", "value"]].rename(
        columns={"value": "affinity_nm"}
    )

    out = presentation.merge(affinity, on=["peptide", "allele"], how="outer")
    return out[[c for c in _SCORE_COLUMNS if c in out.columns]].reset_index(drop=True)


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
    if _predictor_key(predictor) == "mhcflurry":
        return _score_mhcflurry(peptides, alleles)

    try:
        from topiary import TopiaryPredictor
    except ImportError as e:
        raise ImportError(
            "Topiary is required for scoring but is not importable by the Python "
            f"running tsarina ({sys.executable}). Install topiary into that "
            "interpreter, or run tsarina with the environment where topiary is "
            "installed (for example: python -m tsarina.cli panel)."
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
