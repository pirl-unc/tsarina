#!/usr/bin/env python
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Regenerate the tiny real-schema hitlist index in tests/fixtures/hitlist_mini.

``tests/test_hitlist_integration.py`` is the only test that runs the real
``tsarina hits`` path against real hitlist artifacts, so the fixture has to
carry the schema a current hitlist build emits.  Until this script existed the
slice was reproducible only by reverse-engineering its peptide list, which is
how it came to sit a few hitlist releases behind (no ``restriction_evidence``
on observations, no ``gene_biotype`` on mappings).

The slice is taken from the local built index for a handful of MAGEA4 peptides.
Columns the installed hitlist writes but the local artifact predates are filled
in with hitlist's own functions -- the same calls its scanner and builder make
-- so the fixture matches a fresh build without needing one.

Usage
-----
    python scripts/regenerate_hitlist_mini_fixture.py [--check]

``--check`` reports what would change and exits non-zero if anything would,
for use when auditing whether the committed fixture still matches hitlist.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

#: The MAGEA4 peptides the fixture is sliced for. Chosen so the slice exercises
#: multi-mapping paralog attribution (MAGEA4 shares peptides with MAGEA8 and
#: the MAGEA2/2B/12 cluster) while staying a few tens of rows.
FIXTURE_PEPTIDES = (
    "AEMLERVIKNY",
    "AESLFREAL",
    "AETSYVKVL",
    "GVYDGREHTV",
    "WVQENYLEY",
    "YEFLWGPRAL",
)

DEFAULT_OUTPUT_DIR = Path("tests/fixtures/hitlist_mini")


def _annotation_columns(restriction: str) -> dict[str, bool | str]:
    """Every MHC identity column hitlist's scanner persists for one restriction."""
    from hitlist.curation import resolve_mhc_annotation

    return resolve_mhc_annotation(str(restriction)).as_record_fields()


def _restriction_evidence(row: pd.Series) -> str:
    """How the row establishes its named restriction (hitlist#415)."""
    from hitlist.curation import restriction_evidence_for_row

    return restriction_evidence_for_row(
        str(row.get("mhc_restriction", "")),
        pmid=row.get("pmid") if pd.notna(row.get("pmid")) else "",
        submission_id=str(row.get("submission_id", "") or ""),
        is_binding_assay=bool(row.get("is_binding_assay", False)),
        is_monoallelic=bool(row.get("is_monoallelic", False)),
        process_type=str(row.get("process_type", "") or ""),
        disease=str(row.get("disease", "") or ""),
        culture_condition=str(row.get("culture_condition", "") or ""),
        source_tissue=str(row.get("source_tissue", "") or ""),
        cell_name=str(row.get("cell_name", "") or ""),
        assay_comments=str(row.get("assay_comments", "") or ""),
        assay_method=str(row.get("assay_method", "") or ""),
        response_measured=str(row.get("response_measured", "") or ""),
        qualitative_measurement=str(row.get("qualitative_measurement", "") or ""),
    )


def build_observations_slice() -> pd.DataFrame:
    from hitlist.observations import observations_path

    frame = pd.read_parquet(
        observations_path(), filters=[("peptide", "in", list(FIXTURE_PEPTIDES))]
    ).reset_index(drop=True)
    if frame.empty:
        raise SystemExit(
            "No observations rows for the fixture peptides. Build the index first:\n"
            "    tsarina build observations"
        )

    # Re-derive the MHC identity block so every column the current scanner
    # persists is present and internally consistent, however old the artifact
    # this slice came from.
    annotations = pd.DataFrame(
        [_annotation_columns(value) for value in frame["mhc_restriction"]],
        index=frame.index,
    )
    for column in annotations.columns:
        frame[column] = annotations[column]
    frame["restriction_evidence"] = [_restriction_evidence(row) for _, row in frame.iterrows()]
    return frame


def build_mappings_slice() -> pd.DataFrame:
    from hitlist.mappings import mappings_path

    frame = pd.read_parquet(
        mappings_path(), filters=[("peptide", "in", list(FIXTURE_PEPTIDES))]
    ).reset_index(drop=True)
    if frame.empty:
        raise SystemExit(
            "No peptide_mappings rows for the fixture peptides. Build the sidecar first:\n"
            "    tsarina build observations"
        )
    return frame


def _describe(name: str, existing: Path, fresh: pd.DataFrame) -> list[str]:
    """Column and row differences between the committed fixture and a fresh slice."""
    if not existing.exists():
        return [f"{name}: absent, would be written with {len(fresh)} rows"]
    old = pd.read_parquet(existing)
    changes = []
    added = [c for c in fresh.columns if c not in old.columns]
    removed = [c for c in old.columns if c not in fresh.columns]
    if added:
        changes.append(f"{name}: + {', '.join(sorted(added))}")
    if removed:
        changes.append(f"{name}: - {', '.join(sorted(removed))}")
    if len(old) != len(fresh):
        changes.append(f"{name}: {len(old)} rows -> {len(fresh)} rows")
    return changes


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="Fixture directory."
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Report differences without writing; exit 1 if any exist.",
    )
    args = parser.parse_args(argv)

    observations = build_observations_slice()
    mappings = build_mappings_slice()
    targets = {
        "observations.parquet": observations,
        "peptide_mappings.parquet": mappings,
    }

    changes: list[str] = []
    for filename, frame in targets.items():
        changes.extend(_describe(filename, args.output_dir / filename, frame))

    if args.check:
        for line in changes:
            print(line)
        print("fixture matches the installed hitlist" if not changes else "fixture is stale")
        return 1 if changes else 0

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename, frame in targets.items():
        path = args.output_dir / filename
        frame.to_parquet(path, index=False)
        print(f"wrote {path} ({len(frame)} rows x {len(frame.columns)} cols)")
    for line in changes:
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
