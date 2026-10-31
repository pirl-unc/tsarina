"""Current OncoRef decisions projected to gene members for Tsarina consumers."""

from functools import lru_cache

import pandas as pd
from oncoref.cta_selection import results as CTA_selection


@lru_cache(maxsize=1)
def _selection_gene_table() -> pd.DataFrame:
    rows = []
    for decision in CTA_selection(species="human").to_dict("records"):
        identities = {item["gene_id"]: item for item in decision["gene_identities"]}
        for gene_id in decision["gene_ids"]:
            identity = identities[gene_id]
            rows.append(
                {
                    "Ensembl_Gene_ID": gene_id,
                    "Symbol": identity["symbol"] or gene_id,
                    "selection_synonyms": tuple(identity["synonyms"] or ()),
                    **{
                        "selection_" + key: decision[key]
                        for key in (
                            "entity_id",
                            "entity_kind",
                            "label",
                            "gene_ids",
                            "category",
                            "subtype",
                            "assessment_status",
                            "identity_status",
                            "policy_id",
                            "reference_id",
                            "reasons",
                            "evidence_flags",
                            "source_snapshot_id",
                        )
                    },
                }
            )
    frame = pd.DataFrame(rows)
    if frame.Ensembl_Gene_ID.duplicated().any():
        raise ValueError("OncoRef current selection assigns a gene to multiple entities")
    frame["is_current_cta"] = frame.selection_category.eq("CTA").fillna(False).astype(bool)
    return frame


def selection_gene_table() -> pd.DataFrame:
    """Detached gene rows retaining the current protein decision and provenance."""
    return _selection_gene_table().copy()


def selection_gene_ids(category="CTA") -> set[str]:
    frame = _selection_gene_table()
    if category is not None:
        frame = frame.loc[frame.selection_category.eq(category).fillna(False)]
    return set(frame.Ensembl_Gene_ID)


def selection_gene_names(category="CTA") -> set[str]:
    frame = _selection_gene_table()
    if category is not None:
        frame = frame.loc[frame.selection_category.eq(category).fillna(False)]
    return set(frame.Symbol)


def selection_proteoform_groups() -> dict[str, tuple[str, ...]]:
    """Selected reference protein groups with every member gene symbol."""
    frame = _selection_gene_table()
    groups = {}
    for _, members in frame.loc[frame.is_current_cta].groupby("selection_entity_id", sort=False):
        symbols = tuple(sorted(members.Symbol))
        if len(symbols) > 1:
            groups["/".join(symbols)] = symbols
    return groups


@lru_cache(maxsize=1)
def _identity_aliases():
    symbols, synonyms = {}, {}
    for row in _selection_gene_table().itertuples(index=False):
        symbols[row.Symbol.casefold()] = row.Symbol
        for synonym in row.selection_synonyms:
            synonyms.setdefault(synonym.casefold(), set()).add(row.Symbol)
    return symbols, {
        key: next(iter(values)) for key, values in synonyms.items() if len(values) == 1
    }


def selection_symbol_for_alias(name):
    symbols, synonyms = _identity_aliases()
    token = str(name).strip().casefold()
    return symbols.get(token) or synonyms.get(token)


def selection_provenance(*, ensembl_release: int | None = None) -> dict:
    """Selection receipts, independent of a consumer's sequence reference."""
    frame = _selection_gene_table()
    return {
        "cta_selection_policy_ids": sorted(frame.selection_policy_id.dropna().unique()),
        "cta_selection_reference_ids": sorted(frame.selection_reference_id.dropna().unique()),
        "cta_selection_source_snapshot_ids": sorted(
            frame.selection_source_snapshot_id.dropna().unique()
        ),
        "background_policy": "tsarina.current-selection-other-genes.v1",
        **({"ensembl_release": ensembl_release} if ensembl_release is not None else {}),
    }
