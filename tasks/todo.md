# CLI Expansion + Hitlist API Alignment

## Goal

Tsarina's role: **fusion of (a) protein specificity, (b) MS evidence from
hitlist, (c) MHC predictions from topiary**. Deliverables:

1. Fix API mismatches with current hitlist (v1.4.5). Rip out tsarina's
   `hla_only` startswith hack — route everything through hitlist's
   mhcgnomes-based species filter.
2. Make tsarina resolve IEDB/CEDAR paths from hitlist's registry by default.
3. Swap `tsarina.scoring.score_presentation` (direct MHCflurry) to use
   `topiary.TopiaryPredictor` + `mhctools`. MHCflurry becomes one backend
   among several (NetMHCpan, etc.).
4. Add `tsarina personalize` CLI matching `perseus.personalize()`.
5. Add `tsarina hits` CLI: find all public MS hits for a gene/protein, with
   filters for allele / species / serotype / resolution.

## Non-goals

- Rewriting tsarina's own `scan_public_ms` into a thin wrapper over
  `hitlist.scanner.scan`. That's in `curation_tissue_logic_report.md` as a
  separate cleanup; doing it here bundles too much. We will instead make the
  existing wrapper compatible with the new hitlist API and leave consolidation
  for later.
- Fetching IEDB/CEDAR automatically. They remain manual-register (ToS). We just
  surface a clear error pointing at `tsarina data register iedb …`.
- Removing the NetMHCpan `subprocess` path in `scoring.py`. Topiary's NetMHCpan
  backend already wraps it properly; we'll delete `score_netmhcpan` in a later
  pass once callers are migrated.

---

## Phase 1 — Fix breakage + kill `hla_only`

- [x] `tsarina/export.py:71` — `scan(..., hla_only=True)` raises `TypeError`
      against current hitlist. Drop `hla_only`; use
      `mhc_species="Homo sapiens"`.
- [x] Post-scan: filter `hits = hits[~hits.is_binding_assay]` so MS aggregates
      don't silently pick up binding-assay rows (hitlist `1529f66` behavior
      change). Applied in both `export.py` and `evidence.py`.
- [x] `tsarina/iedb.py::scan_public_ms` — removed both `hla_only` *and*
      `human_only`. Single `mhc_species` kwarg (default `"Homo sapiens"`,
      `None` disables) now drives species filtering via
      `hitlist.curation.classify_mhc_species` + `normalize_species`.
- [x] Updated every caller: `perseus.py`, `targets.py`, `viral.py`,
      `mutations.py`, `negatives.py`, `evidence.py`, `export.py`.
- [x] Grep-clean: no `hla_only` or `human_only` in `tsarina/` (test-side
      `startswith("HLA-")` in test_regions/test_alleles is unrelated — it
      validates tsarina's own allele panel definitions).

## Phase 2 — Hitlist as default data source

- [x] New helper `tsarina/datasources.py` with `resolve_iedb_path`,
      `resolve_cedar_path`, `resolve_dataset_paths`, and a
      `DatasetNotRegisteredError` (subclass of `FileNotFoundError`).
- [x] IEDB errors carry the actionable hint
      ``tsarina data register iedb /path/to/mhc_ligand_full.csv``.
      CEDAR silently returns `None` when unregistered.
- [x] Wired into `perseus.personalize`, `export.build_ms_support_maps`,
      `evidence._compute_ms_restriction`. `scan_public_ms` itself stays
      dumb (takes explicit paths); resolution happens at the boundary.
- [x] Added `skip_ms_evidence: bool = False` kwarg to
      `perseus.personalize` to preserve the "don't touch IEDB" path.

## Phase 2b — Topiary for MHC predictions

- [x] Rewrote `tsarina/scoring.py::score_presentation` on top of
      `topiary.TopiaryPredictor` + `mhctools`. `_pivot_topiary` turns
      topiary's long `(peptide, allele, kind)` output into tsarina's wide
      `(peptide, allele, presentation_score, presentation_percentile,
      affinity_nm)` format.
- [x] New `predictor` kwarg — `"mhcflurry"` (default), `"netmhcpan"`,
      `"netmhcpan_el"`.
- [x] `score_affinity` reduced to a thin wrapper over `score_presentation`.
- [x] `score_netmhcpan` kept with a `DeprecationWarning` pointing at
      `score_presentation(predictor="netmhcpan")`.
- [x] `perseus.personalize` now accepts `predictor` and forwards it.

## Phase 3 — CLI `tsarina personalize`

- [ ] New module `tsarina/cli_personalize.py` (argparse subcommand factory +
      handler) to keep `cli.py` readable.
- [ ] Arg surface:
  - `--hla A,B,C,…`                          (required)
  - `--cta GENE=TPM,GENE=TPM`                (repeatable or comma-joined)
  - `--mutations "KRAS G12D,TP53 R175H"`
  - `--viruses hpv16,ebv`
  - `--lengths 8,9,10,11`                    (default 8–11)
  - `--ensembl-release 112`
  - `--mhc-class I|II`                       (default I)
  - `--min-cta-tpm 1.0`
  - `--no-score` (disables topiary scoring; on by default if installed)
  - `--predictor mhcflurry|netmhcpan|netmhcpan_el` (default `mhcflurry`)
  - `--iedb PATH` / `--cedar PATH`           (override registry)
  - `--skip-ms-evidence`                     (don't touch IEDB/CEDAR at all)
  - `-o / --output PATH`                     (CSV; stdout if omitted)
- [ ] Wire `_build_personalize_parser` and `_handle_personalize` into
      `cli.py::main`.

## Phase 4 — CLI `tsarina hits`

A protein → observed-peptide query. Takes a gene symbol (or UniProt ID),
enumerates its k-mers from the Ensembl proteome via
`hitlist.proteome.ProteomeIndex.from_ensembl`, scans IEDB/CEDAR for those
peptides, filters, and aggregates.

- [ ] New module `tsarina/cli_hits.py`.
- [ ] Arg surface:
  - `--gene PRAME` (or `--uniprot P08819`) — one required.
  - `--allele HLA-A*24:02,HLA-A*02:01`       filter on `mhc_restriction`
  - `--species "Homo sapiens"`               pass-through to `scan(mhc_species=…)`
  - `--serotype A2,A24`                      post-scan filter via
    `hitlist.curation.allele_to_serotype`
  - `--min-resolution four_digit|two_digit|serological|class_only`
  - `--mhc-class I|II`
  - `--lengths 8,9,10,11`
  - `--ensembl-release 112`
  - `--predict`                              also run topiary MHC predictions
    for each (peptide, allele) pair; joins `presentation_score`,
    `presentation_percentile`, `affinity_nm` onto output
  - `--predictor mhcflurry|netmhcpan|netmhcpan_el`
  - `--include-binding-assays`               default off
  - `--format peptides|pmhc|raw`             default `pmhc`
    - `peptides`: one row per peptide (via `aggregate_per_peptide`)
    - `pmhc`:     one row per (peptide, allele) (via `aggregate_per_pmhc`)
    - `raw`:      raw scan rows + gene columns from ProteomeIndex
  - `--iedb PATH` / `--cedar PATH`
  - `-o / --output PATH`
- [ ] Gene→peptide resolution strategy:
  1. Build `ProteomeIndex.from_ensembl(release, lengths)` (cached per release).
  2. If `--gene`: scan `idx.protein_meta` for `gene_name == gene`, collect
     `protein_id`s.
  3. If `--uniprot`: direct lookup. (Ensembl releases don't key by UniProt
     accession — if absent, report a clean error pointing the user at
     `--gene`.)
  4. Enumerate all indexed k-mers whose `(protein_id, position)` appears under
     the target protein. Dedupe.
  5. Scan those peptides through hitlist.

## Phase 5 — Verification

- [x] `tests/test_cli_personalize.py` — `--help`, `--skip-ms-evidence
      --no-score` run, and required-arg error all covered.
- [x] `tests/test_cli_hits.py` — `--help`, mutually-exclusive
      `--gene`/`--uniprot`, and "one required" error all covered.
- [x] `tests/test_datasources.py` — monkeypatched `get_path` tests for
      register-hint error, silent-None CEDAR behavior, and
      `require_iedb=False` escape hatch.
- [x] `./format.sh` — 1 file reformatted.
- [x] `./lint.sh` — all checks pass.
- [x] `./test.sh` — 127 passed, 1 skipped (unchanged from baseline).

---

## Review

- **Breakage fixed.** `tsarina/export.py:71` no longer crashes hitlist;
  every scan call passes a supported kwarg set.
- **One species filter.** Removed tsarina's `hla_only` prefix hack and its
  `human_only` host-only check. `mhc_species` (default `"Homo sapiens"`,
  `None` disables) now drives species filtering through
  `hitlist.curation.classify_mhc_species` — single source of truth.
- **Binding-assay drift corrected.** Every scan that backs an MS aggregate
  now filters `~is_binding_assay`. This tightens the MS-only semantics
  `evidence.py` and `export.py` had previously lost.
- **Hitlist is the data manager.** New `datasources.resolve_*` functions are
  the sole place that touches the registry.  CLI, personalize, export, and
  evidence all go through it.  IEDB error message tells the user how to
  register.
- **Topiary swap.** `scoring.score_presentation` wraps
  `topiary.TopiaryPredictor` + `mhctools`, with a `--predictor` switch
  (`mhcflurry` / `netmhcpan` / `netmhcpan_el`). `score_netmhcpan`
  subprocess shim still works but is deprecated.
- **Two new CLIs.**
  - `tsarina personalize` mirrors `perseus.personalize()` kwargs and writes
    a CSV of ranked targets.
  - `tsarina hits` enumerates a gene's k-mers from Ensembl, scans IEDB/CEDAR,
    and filters by allele / species / serotype / resolution / class. Three
    output modes: `peptides` (per-peptide aggregate), `pmhc` (per pMHC),
    `raw` (scan rows plus gene context). `--predict` runs topiary scoring.

## Out-of-scope / follow-ups

- `tsarina/iedb.py::scan_public_ms` still re-implements hitlist's scanner
  logic. Remove this duplication in a later PR (tracked at
  `tasks/curation_tissue_logic_report.md:40`).
- `--uniprot` lookup is best-effort; proper UniProt→Ensembl ID resolution
  needs a mapping file (future).
- Remove `scoring.score_netmhcpan` subprocess shim once topiary-backed
  `score_presentation(predictor="netmhcpan")` proves out.
- hitlist issue #43: `observations.parquet` should always carry
  gene/protein annotations, not collapse multi-mapping to one "best" source.
  Once fixed upstream, `tsarina hits` can swap ProteomeIndex enumeration
  for a simple `load_observations(gene_name=...)` pushdown.

---

## Round 2 — Observations index for fast MS queries

Goal: replace slow CSV rescans with the prebuilt `observations.parquet`
via `hitlist.observations.load_observations`.  Live smoke: 3 peptides
across 90 IEDB rows in 2.6 s end-to-end (vs. multi-minute CSV scan).

- [x] `tsarina/perseus.py` → `tsarina/personalize.py`.  Updated imports in
      `cli_personalize.py`, README, docs/index.md.
- [x] New `tsarina/indexing.py` with `ensure_index_built(force, verbose)`
      and `load_ms_evidence(peptides, mhc_class, mhc_species, ...)`. The
      latter is the one-stop helper: auto-builds, pushdown-filters, and
      applies the peptide / binding-assay filters in memory.
- [x] `tsarina data build [--force]` CLI subcommand.
- [x] Refactored call sites to default to the cached index:
      `personalize.personalize`, `export.build_ms_support_maps`,
      `evidence._compute_ms_restriction`, `cli_hits.handle`.
      Each still honors explicit `--iedb` / `--cedar` overrides by
      falling back to the raw `scan` path.
- [x] `min_allele_resolution` filtering in `cli_hits` moved client-side
      (via `hitlist.curation.classify_allele_resolution` +
      `allele_resolution_rank`) since `load_observations` has no pushdown
      for it.
- [x] `tests/test_indexing.py` — 5 new tests (build skip/trigger/force,
      peptide filter, binding-assay drop). All 132 tests + 1 skip green.
- [x] `./format.sh && ./lint.sh && ./test.sh` all pass.
