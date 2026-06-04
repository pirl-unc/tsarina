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

"""Backfill the ``Aliases`` column of the bundled CTA evidence table.

Per tsarina#77 the column was ~96% empty, so lookups by common antigen names
(``NY-ESO-1`` -> CTAG1B, ``ESO1``, MAGE/GAGE historical symbols) silently
missed. This populates ``Aliases`` for every gene from NCBI ``gene_info``
synonyms (matched by Ensembl gene ID, the authoritative key), merged with any
existing curated aliases and a small supplement of well-known colloquial names
that NCBI does not carry verbatim (e.g. ``NY-ESO-1``).

Source (reproducible):

    curl -L -o human_gene_info.gz \\
      https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz

Run ``python scripts/backfill_aliases.py --check`` to preview coverage without
writing.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
CSV_PATH = REPO_ROOT / "tsarina" / "data" / "cancer-testis-antigens.csv"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cta_sources import ncbi_gene_info_path  # noqa: E402

#: Well-known colloquial CTA names not present verbatim in NCBI synonyms,
#: keyed by official symbol.  Keep small and high-confidence.
_CURATED_COLLOQUIAL: dict[str, list[str]] = {
    "CTAG1B": ["NY-ESO-1"],
    "CTAG2": ["LAGE-1", "CT6.2"],
}


def _load_ncbi_synonyms(gene_info_path: Path) -> dict[str, list[str]]:
    """Map (unversioned) Ensembl gene ID -> NCBI synonym list."""
    ensg_to_syn: dict[str, list[str]] = {}
    with gzip.open(gene_info_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        i_syn = header.index("Synonyms")
        i_xref = header.index("dbXrefs")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            raw = parts[i_syn]
            synonyms = [] if raw == "-" else [s for s in raw.split("|") if s and s != "-"]
            for xref in parts[i_xref].split("|"):
                if xref.startswith("Ensembl:"):
                    ensg_to_syn[xref.split(":", 1)[1]] = synonyms
    return ensg_to_syn


def _split_existing(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [a.strip() for a in str(value).split(";") if a.strip() and a.strip().lower() != "nan"]


def build_aliases(df: pd.DataFrame, ensg_to_syn: dict[str, list[str]]) -> pd.Series:
    out = []
    for _, row in df.iterrows():
        symbol = str(row["Symbol"]).strip()
        ensg = str(row["Ensembl_Gene_ID"]).split(".")[0]
        merged: list[str] = []
        seen_norm: set[str] = {symbol.upper()}
        for alias in (
            _split_existing(row.get("Aliases"))
            + ensg_to_syn.get(ensg, [])
            + _CURATED_COLLOQUIAL.get(symbol, [])
        ):
            norm = alias.upper()
            if norm in seen_norm:
                continue
            seen_norm.add(norm)
            merged.append(alias)
        out.append(";".join(sorted(merged)) if merged else "")
    return pd.Series(out, index=df.index)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gene-info",
        default=None,
        help="Local NCBI gene_info.gz path; default downloads + caches from NCBI.",
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    df = pd.read_csv(CSV_PATH)
    ensg_to_syn = _load_ncbi_synonyms(ncbi_gene_info_path(args.gene_info))
    new_aliases = build_aliases(df, ensg_to_syn)

    nonempty = int((new_aliases.str.len() > 0).sum())
    print(f"Genes with >=1 alias: {nonempty}/{len(df)} "
          f"(was {int(df['Aliases'].notna().sum())})")
    for sym in ("CTAG1B", "MAGEA4", "PRAME", "XAGE1B"):
        idx = df.index[df["Symbol"] == sym]
        if len(idx):
            print(f"  {sym}: {new_aliases.loc[idx[0]]}")

    if args.check:
        print("\n--check: no files written.")
        return

    df["Aliases"] = new_aliases
    df.to_csv(CSV_PATH, index=False)
    print(f"\nWrote {len(df)} rows to {CSV_PATH}")


if __name__ == "__main__":
    main()
