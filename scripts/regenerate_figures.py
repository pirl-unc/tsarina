#!/usr/bin/env python
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Regenerate every tsarina figure family into one timestamped run.

Writes ``figures/run_<YYYYMMDD-HHMMSS>/{cta,neoantigen-vaccine}/*.png``,
bundles every PNG into ``figures/run_.../tsarina-all-figures.pdf``, and repoints
the top-level ``figures/latest`` symlink at the new run -- the same
convention independently adopted for oncoref's figure output, so a
combined multi-repo ``figures/`` tree reads as one system.

Usage::

    python scripts/regenerate_figures.py
    python scripts/regenerate_figures.py --out-dir figures/my-run
"""

from __future__ import annotations

import argparse
import time
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")


def _resolve_run_dir(explicit: Path | None) -> Path:
    if explicit is not None:
        explicit.mkdir(parents=True, exist_ok=True)
        return explicit
    root = Path("figures")
    stamp = time.strftime("run_%Y%m%d-%H%M%S")
    run_dir = root / stamp
    suffix = 1
    while run_dir.exists():
        suffix += 1
        run_dir = root / f"{stamp}-{suffix}"
    run_dir.mkdir(parents=True)
    return run_dir


def _write_all_figures_pdf(run_dir: Path, written: list[Path]) -> Path | None:
    from matplotlib.backends.backend_pdf import PdfPages

    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt

    if not written:
        return None
    pdf_path = run_dir / "tsarina-all-figures.pdf"
    with PdfPages(pdf_path) as pages:
        for path in written:
            img = mpimg.imread(path)
            h, w = img.shape[0], img.shape[1]
            fig = plt.figure(figsize=(w / 150, h / 150), dpi=150)
            ax = fig.add_axes((0, 0, 1, 1))
            ax.imshow(img)
            ax.axis("off")
            pages.savefig(fig)
            plt.close(fig)
    return pdf_path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None, help="base output dir (default: figures/run_<timestamp>)")
    ap.add_argument(
        "--vaccine-peptide-table",
        type=Path,
        default=Path("tests/data/osteosarc/vaccine_overlap_summary.tsv"),
        help="wide-format vaccine-overlap-summary TSV for the neoantigen family",
    )
    args = ap.parse_args()

    run_dir = _resolve_run_dir(args.out_dir)
    print(f"Regenerating tsarina figures into {run_dir}")

    written: list[Path] = []
    skipped: list[tuple[str, str]] = []

    try:
        from tsarina.ms_evidence_plots import render_all as render_cta

        cta_dir = run_dir / "cta"
        result = render_cta(out_dir=cta_dir)
        for name, path in result.items():
            written.append(path)
            print(f"  ok    cta/{path.name}")
    except Exception as e:
        skipped.append(("cta", f"{type(e).__name__}: {e}"))
        print(f"  SKIP  cta  ({type(e).__name__}: {e})")
        traceback.print_exc()

    try:
        from tsarina.neoantigen_evidence_plots import render_all as render_neoantigen

        neoantigen_dir = run_dir / "neoantigen-vaccine"
        result = render_neoantigen(args.vaccine_peptide_table, out_dir=neoantigen_dir)
        for name, path in result.items():
            written.append(path)
            print(f"  ok    neoantigen-vaccine/{path.name}")
    except Exception as e:
        skipped.append(("neoantigen-vaccine", f"{type(e).__name__}: {e}"))
        print(f"  SKIP  neoantigen-vaccine  ({type(e).__name__}: {e})")
        traceback.print_exc()

    pdf = _write_all_figures_pdf(run_dir, written)
    if pdf is not None:
        print(f"  ok    {pdf.name}")

    if args.out_dir is None:
        latest = Path("figures") / "latest"
        if latest.is_symlink() or latest.exists():
            latest.unlink()
        latest.symlink_to(run_dir.name)
        print(f"figures/latest -> {run_dir.name}")

    print(f"\n{len(written)} figures written, {len(skipped)} family(ies) skipped.")
    for family, reason in skipped:
        print(f"  {family}: {reason}")
    return 1 if skipped and not written else 0


if __name__ == "__main__":
    raise SystemExit(main())
