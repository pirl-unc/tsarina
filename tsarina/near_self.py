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

"""OncoRef 1.52+ CTA safety caveats and 1-aa Hamming near-self screening.

Two complementary layers of TCR cross-reactivity / target-safety visibility:

1. :func:`annotate_cta_safety_caveats` — instant (<1 ms) per-row annotation from
   OncoRef's target-selection metadata: ``nuclear_rna_caution`` (intronic snRNA-seq
   or bidirectional-promoter spillover), epitope-scoped clinical organ clearance
   (distinguishing cleared epitopes such as ``CTAG2`` ``SLLMWITQC`` from unverified
   epitopes on the same locus), and high 1-mm near-self overlap
   (``unsafe_near_self_1mm_pct >= 40%``).
2. :func:`screen_near_self_1mm` / :func:`annotate_near_self_1mm` — opt-in
   per-peptide 1-aa Hamming-distance (substitution-only, no indels) screen
   against OncoRef's ``Unsafe`` (or ``Unsafe|Other``) proteome via
   :class:`oncoref.peptide_subtraction.PeptideIndex`.
"""

from __future__ import annotations

import json
import math
from collections.abc import Iterable
from functools import lru_cache

import pandas as pd
from oncoref.peptide_subtraction import (
    AMINO_ACIDS,
    NEAR_SELF_CONTRACT_ID,
    MatchPolicy,
    PeptideIndex,
)

#: Gene-level threshold (%) above which ``cta_safety_note`` flags high 1-mm near-self risk.
HIGH_NEAR_SELF_1MM_PCT_THRESHOLD: float = 40.0

#: Supported OncoRef exclusion backgrounds for 1-mm near-self screening.
VALID_NEAR_SELF_BACKGROUNDS: tuple[str, ...] = ("Unsafe", "Unsafe|Other")

#: Exact 1-aa Hamming-distance MatchPolicy (1 substitution allowed, indels disallowed).
HAMMING_1MM_POLICY = MatchPolicy(
    max_distance=1.0,
    insertion_cost=2.0,
    deletion_cost=2.0,
)


def _parse_json_list(value: object) -> tuple[str, ...]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ()
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped or stripped.lower() == "nan":
            return ()
        if stripped.startswith("["):
            try:
                parsed = json.loads(stripped)
            except json.JSONDecodeError:
                return (stripped,)
            return tuple(str(item).strip() for item in parsed if str(item).strip())
        return (stripped,)
    if isinstance(value, Iterable):
        return tuple(str(item).strip() for item in value if str(item).strip())
    return (str(value).strip(),)


def _clean_str(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def _clean_float(value: object) -> float | None:
    if value is None:
        return None
    try:
        num = float(value)
    except (TypeError, ValueError):
        return None
    return num if math.isfinite(num) else None


@lru_cache(maxsize=8)
def _cta_safety_lookup(species_key: str) -> dict[str, dict[str, object]]:
    """Map CTA display labels and member symbols to OncoRef safety metadata."""
    from .cta_selection import CTA_selection

    sel = CTA_selection(species=species_key)
    if sel.empty:
        return {}

    proteoform_symbol_fn = None
    if species_key in ("human", "homo_sapiens", "homo sapiens", "9606"):
        try:
            from oncoref.proteoforms import proteoform_symbol as _ps

            proteoform_symbol_fn = _ps
        except Exception:
            proteoform_symbol_fn = None

    lookup: dict[str, dict[str, object]] = {}
    for row in sel.itertuples(index=False):
        label = _clean_str(getattr(row, "label", ""))
        if not label:
            continue
        info: dict[str, object] = {
            "label": label,
            "rare_cell_status": _clean_str(getattr(row, "rare_cell_status", "")),
            "nuclear_rna_caution": _clean_str(getattr(row, "nuclear_rna_caution", "")),
            "clinical_clearance_scope": _clean_str(getattr(row, "clinical_clearance_scope", "")),
            "clinical_cleared_organs": _parse_json_list(
                getattr(row, "clinical_cleared_organs", None)
            ),
            "clinical_cleared_epitopes": frozenset(
                e.upper() for e in _parse_json_list(getattr(row, "clinical_cleared_epitopes", None))
            ),
            "unsafe_near_self_1mm_pct": _clean_float(
                getattr(row, "unsafe_near_self_1mm_pct", None)
            ),
        }
        keys = {label}
        if "/" in label:
            keys.update(part.strip() for part in label.split("/") if part.strip())
        if proteoform_symbol_fn is not None:
            try:
                preferred = _clean_str(proteoform_symbol_fn(label))
                if preferred:
                    keys.add(preferred)
            except Exception:
                pass
        identities = getattr(row, "gene_identities", None)
        if isinstance(identities, (tuple, list)):
            for item in identities:
                if isinstance(item, dict):
                    sym = _clean_str(item.get("symbol", ""))
                    if sym:
                        keys.add(sym)
        for key in keys:
            # Prefer CTA rows over non-CTA rows if an alias collides.
            existing = lookup.get(key)
            if existing is None or getattr(row, "category", None) == "CTA":
                lookup[key] = info
    return lookup


def _resolve_safety_info(
    source: str,
    lookup: dict[str, dict[str, object]],
) -> dict[str, object] | None:
    if not source:
        return None
    direct = lookup.get(source)
    if direct is not None:
        return direct
    if "/" in source:
        for part in source.split("/"):
            hit = lookup.get(part.strip())
            if hit is not None:
                return hit
    return None


def _build_safety_note(
    info: dict[str, object],
    peptide: str,
) -> tuple[str, bool, float | None]:
    cleared_epitopes: frozenset[str] = info["clinical_cleared_epitopes"]  # type: ignore[assignment]
    pep_upper = peptide.strip().upper()
    is_cleared = bool(pep_upper and pep_upper in cleared_epitopes)
    near_self_pct: float | None = info["unsafe_near_self_1mm_pct"]  # type: ignore[assignment]

    notes: list[str] = []
    nuclear_caution = str(info.get("nuclear_rna_caution") or "")
    if nuclear_caution:
        notes.append(f"nuclear RNA caution ({nuclear_caution})")

    clearance_scope = str(info.get("clinical_clearance_scope") or "")
    rare_cell = str(info.get("rare_cell_status") or "")
    if (
        clearance_scope == "epitope"
        and rare_cell == "unresolved"
        and cleared_epitopes
        and not is_cleared
    ):
        epitope_list = ", ".join(sorted(cleared_epitopes))
        notes.append(f"clinical clearance scoped to {epitope_list} (unverified epitope)")

    if near_self_pct is not None and near_self_pct >= HIGH_NEAR_SELF_1MM_PCT_THRESHOLD:
        notes.append(
            f"high 1-mm near-self overlap ({near_self_pct:.1f}% of exact-unique 9-mers "
            "within 1 aa of Unsafe proteome)"
        )

    return "; ".join(notes), is_cleared, near_self_pct


def annotate_cta_safety_caveats(
    df: pd.DataFrame,
    *,
    source_col: str = "source",
    peptide_col: str = "peptide",
    category_col: str | None = "category",
    species: str = "human",
) -> pd.DataFrame:
    """Attach OncoRef 1.52+ gene/epitope CTA safety caveats onto a target table.

    Adds three columns without dropping rows:

    - ``cta_safety_note`` (str): concise caveat when ``nuclear_rna_caution`` is
      set, when clinical clearance is epitope-scoped on an unresolved locus and
      the row's peptide is not among the cleared epitopes (e.g., non-``SLLMWITQC``
      ``CTAG2`` peptides), or when ``unsafe_near_self_1mm_pct >= 40%``. Empty
      string otherwise.
    - ``clinical_cleared_epitope`` (bool): ``True`` when the peptide exactly
      matches an OncoRef clinically cleared epitope for that CTA (e.g.,
      ``CTAG2`` ``SLLMWITQC``, ``MAGEA4`` ``GVYDGREHTV``, ``PRAME`` ``SLLQHLIGL``).
    - ``unsafe_near_self_1mm_pct`` (float / NA): OncoRef's percentage of
      exact-unique 9-mers within 1 aa of the ``Unsafe`` proteome.
    """
    out = df.copy()
    if out.empty:
        out["unsafe_near_self_1mm_pct"] = pd.Series(dtype="float64")
        out["clinical_cleared_epitope"] = pd.Series(dtype="bool")
        out["cta_safety_note"] = pd.Series(dtype="object")
        return out

    species_key = str(species).strip().lower()
    lookup = _cta_safety_lookup(species_key)

    notes: list[str] = []
    cleared_flags: list[bool] = []
    pct_values: list[float | None] = []

    has_cat = category_col is not None and category_col in out.columns
    sources = out[source_col].fillna("").astype(str) if source_col in out.columns else None
    peptides = out[peptide_col].fillna("").astype(str) if peptide_col in out.columns else None
    categories = out[category_col].fillna("").astype(str) if has_cat else None

    for idx in range(len(out)):
        if categories is not None and categories.iloc[idx] not in {"cta", "cta_flagged"}:
            notes.append("")
            cleared_flags.append(False)
            pct_values.append(None)
            continue
        source = sources.iloc[idx].strip() if sources is not None else ""
        peptide = peptides.iloc[idx].strip() if peptides is not None else ""
        info = _resolve_safety_info(source, lookup)
        if info is None:
            notes.append("")
            cleared_flags.append(False)
            pct_values.append(None)
            continue
        note, is_cleared, pct = _build_safety_note(info, peptide)
        notes.append(note)
        cleared_flags.append(is_cleared)
        pct_values.append(pct)

    out["unsafe_near_self_1mm_pct"] = pd.Series(pct_values, index=out.index, dtype="float64")
    out["clinical_cleared_epitope"] = pd.Series(cleared_flags, index=out.index, dtype="bool")
    out["cta_safety_note"] = pd.Series(notes, index=out.index, dtype="object")
    return out


def _normalize_background_proteins(
    background_proteins: pd.DataFrame | dict[str, str] | Iterable[str],
) -> pd.DataFrame:
    """Normalize caller-supplied background proteins into a deduplicated PeptideIndex frame."""
    seq_to_symbols: dict[str, set[str]] = {}
    if isinstance(background_proteins, pd.DataFrame):
        if "sequence" not in background_proteins.columns:
            raise ValueError("background_proteins DataFrame must contain a 'sequence' column")
        sym_col = (
            "symbol"
            if "symbol" in background_proteins.columns
            else "gene_name"
            if "gene_name" in background_proteins.columns
            else "gene_id"
            if "gene_id" in background_proteins.columns
            else None
        )
        for idx, row in enumerate(background_proteins.itertuples(index=False)):
            seq = _clean_str(getattr(row, "sequence", "")).upper()
            if not seq:
                continue
            sym = _clean_str(getattr(row, sym_col, "")) if sym_col else f"bg_{idx}"
            seq_to_symbols.setdefault(seq, set()).add(sym or f"bg_{idx}")
    elif isinstance(background_proteins, dict):
        for label, seq_raw in background_proteins.items():
            seq = _clean_str(seq_raw).upper()
            if seq:
                seq_to_symbols.setdefault(seq, set()).add(str(label).strip() or "bg")
    else:
        for idx, seq_raw in enumerate(background_proteins):
            seq = _clean_str(seq_raw).upper()
            if seq:
                seq_to_symbols.setdefault(seq, set()).add(f"bg_{idx}")

    if not seq_to_symbols:
        return pd.DataFrame(columns=["symbol", "sequence"])
    return pd.DataFrame(
        [
            {"symbol": ";".join(sorted(syms)), "sequence": seq}
            for seq, syms in seq_to_symbols.items()
        ]
    )


def _background_gene_ids(species_key: str, background: str) -> set[str]:
    from .cta_selection import CTA_selection

    categories = ("Unsafe",) if background == "Unsafe" else ("Unsafe", "Other")
    sel = CTA_selection(category=categories, species=species_key)
    gene_ids: set[str] = set()
    for members in sel.get("gene_ids", ()):
        if isinstance(members, (tuple, list)):
            gene_ids.update(str(gid).split(".")[0] for gid in members if gid)
    return gene_ids


@lru_cache(maxsize=4)
def _cached_proteome_1mm_index(
    species_key: str,
    ensembl_release: int,
    background: str,
    lengths: tuple[int, ...],
) -> PeptideIndex:
    """Build and cache a 1-aa Hamming PeptideIndex over OncoRef Unsafe[/Other] Ensembl proteins."""
    from pyensembl import EnsemblRelease

    from .species import resolve_species_context

    ctx = resolve_species_context(species_key, ensembl_release=ensembl_release)
    bg_ids = _background_gene_ids(ctx.host_species, background)
    if not bg_ids:
        return PeptideIndex(
            pd.DataFrame(columns=["symbol", "sequence"]),
            policy=HAMMING_1MM_POLICY,
            lengths=lengths,
        )

    ensembl = (
        EnsemblRelease(ensembl_release)
        if ctx.is_human_host
        else EnsemblRelease(ensembl_release, species=ctx.pyensembl_species)
    )
    fasta_by_id = {
        str(k).split(".")[0]: v for k, v in ensembl.protein_sequences.fasta_dictionary.items()
    }
    cur = ensembl.db.connection.execute(
        "SELECT DISTINCT gene_id, gene_name, protein_id FROM CDS "
        "WHERE protein_id IS NOT NULL AND protein_id != ''"
    )
    seq_to_symbols: dict[str, set[str]] = {}
    for gene_id, gene_name, protein_id in cur.fetchall():
        gid = str(gene_id).split(".")[0]
        if gid not in bg_ids:
            continue
        seq = fasta_by_id.get(str(protein_id).split(".")[0])
        if not seq:
            continue
        sym = str(gene_name).strip() if gene_name else gid
        seq_to_symbols.setdefault(seq, set()).add(sym or gid)

    records = pd.DataFrame(
        [
            {"symbol": ";".join(sorted(syms)), "sequence": seq}
            for seq, syms in seq_to_symbols.items()
        ],
        columns=["symbol", "sequence"],
    )
    return PeptideIndex(records, policy=HAMMING_1MM_POLICY, lengths=lengths)


def screen_near_self_1mm(
    peptides: Iterable[str],
    *,
    ensembl_release: int = 112,
    species: str = "human",
    background: str = "Unsafe",
    background_proteins: pd.DataFrame | dict[str, str] | Iterable[str] | None = None,
) -> pd.DataFrame:
    """Screen candidate peptides for 1-aa Hamming-distance matches in the Unsafe proteome.

    Uses :class:`oncoref.peptide_subtraction.PeptideIndex` with exact 1-aa
    substitution distance (no indels; exact 0-mm hits are excluded from the
    1-mm count since exact self-matches are screened separately).

    Returns a DataFrame with one row per unique input peptide:
    ``peptide``, ``near_self_1mm_count``, ``near_self_1mm_genes``,
    ``near_self_1mm_peptides``.
    """
    if background not in VALID_NEAR_SELF_BACKGROUNDS:
        raise ValueError(
            f"near_self_background must be one of {VALID_NEAR_SELF_BACKGROUNDS}, got {background!r}"
        )

    unique_peps = list(dict.fromkeys(str(p).strip().upper() for p in peptides if str(p).strip()))
    columns = [
        "peptide",
        "near_self_1mm_count",
        "near_self_1mm_genes",
        "near_self_1mm_peptides",
    ]
    if not unique_peps:
        empty = pd.DataFrame(columns=columns)
        empty.attrs.update(
            near_self_contract_id=NEAR_SELF_CONTRACT_ID,
            near_self_background=background,
            near_self_species=species,
            near_self_ensembl_release=ensembl_release,
        )
        return empty

    lengths = tuple(sorted({len(p) for p in unique_peps}))
    if background_proteins is not None:
        bg_df = _normalize_background_proteins(background_proteins)
        index = (
            PeptideIndex(bg_df, policy=HAMMING_1MM_POLICY, lengths=lengths)
            if not bg_df.empty
            else None
        )
    else:
        index = _cached_proteome_1mm_index(
            str(species).strip().lower(),
            int(ensembl_release),
            background,
            lengths,
        )

    rows: list[dict[str, object]] = []
    for pep in unique_peps:
        if index is None or not set(pep) <= AMINO_ACIDS:
            rows.append(
                {
                    "peptide": pep,
                    "near_self_1mm_count": 0,
                    "near_self_1mm_genes": "",
                    "near_self_1mm_peptides": "",
                }
            )
            continue
        matches = index.matches(pep)
        if matches.empty:
            rows.append(
                {
                    "peptide": pep,
                    "near_self_1mm_count": 0,
                    "near_self_1mm_genes": "",
                    "near_self_1mm_peptides": "",
                }
            )
            continue
        # Restrict strictly to 1-aa substitutions (distance > 0, same length).
        mm1 = matches.loc[
            (matches["distance"] > 0) & (matches["matched_peptide"].str.len() == len(pep))
        ]
        if mm1.empty:
            rows.append(
                {
                    "peptide": pep,
                    "near_self_1mm_count": 0,
                    "near_self_1mm_genes": "",
                    "near_self_1mm_peptides": "",
                }
            )
            continue
        matched_peps = sorted(set(mm1["matched_peptide"].astype(str)))
        matched_genes: set[str] = set()
        if "symbol" in mm1.columns:
            for cell in mm1["symbol"].dropna().astype(str):
                for part in cell.split(";"):
                    token = part.strip()
                    if token:
                        matched_genes.add(token)
        rows.append(
            {
                "peptide": pep,
                "near_self_1mm_count": len(matched_peps),
                "near_self_1mm_genes": ";".join(sorted(matched_genes)),
                "near_self_1mm_peptides": ";".join(matched_peps),
            }
        )

    out = pd.DataFrame(rows, columns=columns)
    out.attrs.update(
        near_self_contract_id=NEAR_SELF_CONTRACT_ID,
        near_self_background=background,
        near_self_species=species,
        near_self_ensembl_release=ensembl_release,
    )
    return out


def annotate_near_self_1mm(
    df: pd.DataFrame,
    *,
    peptide_col: str = "peptide",
    ensembl_release: int = 112,
    species: str = "human",
    background: str = "Unsafe",
    background_proteins: pd.DataFrame | dict[str, str] | Iterable[str] | None = None,
    max_matches: int | None = None,
) -> pd.DataFrame:
    """Attach 1-aa near-self match counts/genes/peptides onto a target table and optionally gate."""
    if max_matches is not None and max_matches < 0:
        raise ValueError(f"max_matches must be >= 0, got {max_matches!r}")

    out = df.copy()
    if out.empty or peptide_col not in out.columns:
        out["near_self_1mm_count"] = pd.Series(dtype="int64")
        out["near_self_1mm_genes"] = pd.Series(dtype="object")
        out["near_self_1mm_peptides"] = pd.Series(dtype="object")
        return out

    screened = screen_near_self_1mm(
        out[peptide_col].dropna().astype(str).unique(),
        ensembl_release=ensembl_release,
        species=species,
        background=background,
        background_proteins=background_proteins,
    )
    out = out.drop(
        columns=[
            c
            for c in ("near_self_1mm_count", "near_self_1mm_genes", "near_self_1mm_peptides")
            if c in out.columns
        ],
        errors="ignore",
    )
    out["_near_self_pep_key"] = out[peptide_col].astype(str).str.strip().str.upper()
    out = out.merge(
        screened.rename(columns={"peptide": "_near_self_pep_key"}),
        on="_near_self_pep_key",
        how="left",
        validate="m:1",
    ).drop(columns=["_near_self_pep_key"])
    out["near_self_1mm_count"] = (
        pd.to_numeric(out["near_self_1mm_count"], errors="coerce").fillna(0).astype(int)
    )
    out["near_self_1mm_genes"] = out["near_self_1mm_genes"].fillna("").astype(str)
    out["near_self_1mm_peptides"] = out["near_self_1mm_peptides"].fillna("").astype(str)
    if max_matches is not None:
        out = out.loc[out["near_self_1mm_count"] <= int(max_matches)].reset_index(drop=True)
    out.attrs.update(screened.attrs)
    return out
