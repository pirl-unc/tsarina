# PR — Panel CTA Safety And NY-ESO-1 Grouping (2026-04-30)

## Goal

Make `tsarina panel` select a less MAGE-heavy, clinically anchored CTA set by
grouping `CTAG1A`/`CTAG1B` as one `NY-ESO-1` target and filtering automatic CTA
selection against vital-tissue RNA expression unless explicitly allowed.

## Plan

- [x] Add CTA alias/group resolution so `NY-ESO-1` expands to `CTAG1A` and
      `CTAG1B`, then collapses peptide rows back to a single display CTA.
- [x] Add an automatic vital-tissue RNA/MS gate using existing CTA table
      columns, with default allow-list exceptions for `PRAME`, `NY-ESO-1`, and
      `MAGEA4`.
- [x] Expose CLI parameters for the vital-tissue gate and allow-list.
- [x] Document the updated panel defaults and alias handling.
- [x] Add tests for PRAME inclusion, NY-ESO-1 grouping, MAGE alias handling,
      and CLI wiring.
- [x] Bump patch version and run full verification.

## Verification

- [x] `./format.sh`
- [x] `./lint.sh`
- [x] `./test.sh` — 254 passed

# PR — Spanning Alias CLI Smoke Tests (2026-04-30)

## Goal

Close #23 by pinning down the deprecated `tsarina spanning` CLI alias with
subprocess smoke tests equivalent to the visible `tsarina panel` command tests.

## Plan

- [x] Add headline `spanning --help` assertions for the requested flags.
- [x] Add invalid `--panel` and `--format` tests for the alias.
- [x] Bump patch version and run full verification.

## Verification

- [x] `./format.sh`
- [x] `./lint.sh`
- [x] `./test.sh` — 251 passed

# PR — Panel Version And Peptide Progress (2026-04-30)

## Goal

Make `tsarina panel` show the package version and enough progress to explain
slow peptide enumeration. Remove the current startup bottleneck where CTA
exclusivity builds a full human non-CTA k-mer posting index before any useful
panel progress appears.

## Plan

- [x] Print `tsarina vX.Y.Z` at the start of panel progress output.
- [x] Thread panel progress callbacks into peptide generation.
- [x] Limit panel peptide enumeration to the selected CTA list instead of
      generating all CTAs first.
- [x] Replace CTA exclusivity's full non-CTA `proteome_kmer_set()` build with
      a streaming non-CTA protein scan against the much smaller CTA peptide set.
- [x] Add progress messages and optional tqdm bars for CTA gene enumeration and
      non-CTA background scanning.
- [x] Make the Topiary missing-import error name the Python interpreter that
      cannot import it.
- [x] Add tests for version output, peptide-stage progress, and the streaming
      exclusivity filter.
- [x] Bump patch version and run full verification.

## Verification

- [x] `./format.sh`
- [x] `./lint.sh`
- [x] `./test.sh` — 249 passed

# PR — Panel Progress, Table Output, And Summary (2026-04-30)

## Goal

Make `tsarina panel` useful as an interactive default command, not just a CSV
generator: show meaningful progress while long scoring work runs, print a
readable terminal table by default, keep CSV output for automation, and append
coverage-style summary statistics.

## Plan

- [x] Add richer panel progress stages for CTA resolution, peptide enumeration,
      public-MS evidence loading, scoring, candidate selection, and summary.
- [x] Score in chunks when CLI progress bars are enabled so `tqdm` can report
      the slow prediction stage instead of a single quiet blocking call.
- [x] Support a configurable top-N peptide cap per CTA x HLA cell, default 3,
      ranked by MS source count first, then MS hit count, then prediction.
- [x] Preserve current wide/long CSV modes while adding a default text table
      CLI mode with CTA rows ordered by CTA rank and HLA alleles ordered by
      weighted population frequency when built-in frequency data exist.
- [x] Add summary statistics for HLA allele count, CTA count, selected peptide
      count, filled cells, expected population coverage per CTA, and fraction
      of CTAs covered per HLA.
- [x] Update CLI help/docs and bump the package patch version.
- [x] Add focused tests for default CLI table output, progress behavior,
      top-N selection/ranking, weighted ordering, and summary calculation.

## Verification

- [x] `./format.sh`
- [x] `./lint.sh`
- [x] `./test.sh` — 246 passed

# PR — Panel Command And MS Evidence Tiers (2026-04-30)

## Goal

Replace the user-facing `spanning` command with a clearer `panel` command that
builds a CTA x HLA pMHC matrix from sensible defaults, while making panel cell
selection MS-evidence-first and exposing configurable evidence-tier prediction
cutoffs.

## Plan

- [x] Rename the primary CLI command to `panel` with help text:
      "Build a CTA x HLA pMHC matrix for a population HLA panel."
- [x] Keep `spanning` as a deprecated compatibility alias for one release.
- [x] Change panel defaults to top 25 CTAs, `global51_abc_ssa`, lengths
      `8,9,10,11`, MHCflurry, and MS-evidence-first selection.
- [x] Classify candidate peptide-HLA cells into `monoallelic_ms`,
      `sample_allele_ms`, `unrestricted_ms`, and optional `predicted_only`.
- [x] Use default tier cutoffs of 2.0, 1.0, 0.5, and 0.1 percentile,
      respectively, and expose them as library parameters and CLI flags.
- [x] Ensure `predicted_only` is excluded by default and, when enabled, cannot
      displace an MS-supported candidate.
- [x] Add focused tests for defaults, tier classification, configurable
      thresholds, CLI flags, and deprecated alias behavior.

## Verification

- [x] `./format.sh`
- [x] `./lint.sh`
- [x] `./test.sh` — 240 passed

# PR Series — Evidence And CLI Cleanup (2026-04-30)

## Goal

Start a stacked set of small PRs that fixes the audit issues while reducing
duplicate call paths. Keep each PR reviewable, with one version bump per PR.

## Planned Stack

- [x] PR 1: Centralize public-MS loading and HLA normalization helpers; fix
      `target_peptides()` default evidence filtering (#31) and
      `tsarina hits --allele` normalization (#33).
- [x] PR 2: Rework panel matrix metrics to use the centralized helpers, fail
      loudly when scoring is unavailable, and count literal HLA restrictions
      instead of regex strings (#34).
- [x] PR 3: Refresh README/docs around hitlist data location and the default
      human-exclusive viral helper (#32).

## Design Notes

- Add one public-MS loader that chooses the cached hitlist observations path by
  default and the raw scanner path only when explicit IEDB/CEDAR paths are
  supplied.
- Keep binding-assay exclusion in that loader so target, export, evidence, and
  panel paths cannot diverge.
- Add shared MHC restriction normalization/matching helpers and call them from
  CLI filters and panel MS metrics.
- Preserve existing output columns where callers may depend on them; add tests
  for each fixed behavior before opening PRs.

## PR 1 Verification

- [x] `./format.sh`
- [x] `./lint.sh`
- [x] `./test.sh` — 229 passed

## PR 2 Verification

- [x] `./format.sh`
- [x] `./lint.sh`
- [x] `./test.sh` — 235 passed

## PR 3 Verification

- [x] `./format.sh`
- [x] `./lint.sh`
- [x] `./test.sh` — 235 passed

# Audit — Hitlist Changes And Full Tsarina Pass (2026-04-29)

## Goal

Inspect the current hitlist state that tsarina depends on, then read through
tsarina end to end looking for bugs, API mismatches, stale assumptions,
documentation drift, test gaps, and simple improvement opportunities.

## Plan

- [x] Inspect local hitlist checkout, installed hitlist version, and recent
      hitlist API changes relevant to tsarina.
- [x] Review tsarina package metadata, public modules, CLIs, tests, and task
      notes so the audit is grounded in current code rather than prior reports.
- [x] Trace hitlist-facing paths: data registry, observation loading, scanning,
      aggregation, proteome k-mer delegation, and CLI filters.
- [x] Record concrete findings with file/line references and classify them by
      severity.
- [x] Decide whether any findings are clearly safe to fix now; otherwise leave
      an actionable review section.

## Review

### Context Checked

- Local `hitlist` import points at `/Users/iskander/code/hitlist`; package
  attribute reports `1.30.4`, while installed metadata reports `1.30.0`.
- Sibling `hitlist` checkout is on `feat/pmhc-sample-paired` with one
  modified file: `hitlist/pmhc_query.py`.
- Current `tsarina` main includes PR #30 / v1.2.0:
  hitlist `>=1.15.1`, `tsarina hits --lengths` pushdown, and
  class-derived length defaults.

### Findings

1. **High — `target_peptides()` ignores the cached hitlist evidence path.**
   `tsarina/targets.py:177-249` only attaches MS evidence when explicit
   `iedb_path` or `cedar_path` is supplied. With defaults, both
   `require_ms_evidence=True` and `cancer_specific=True` are no-ops.
   When explicit raw CSVs are supplied, binding-assay rows are not dropped
   before MS aggregation. Filed as
   https://github.com/pirl-unc/tsarina/issues/31.

2. **High — panel MS metrics are not reliable.**
   `tsarina/panels.py:113-188` inherits the `target_peptides()` evidence
   no-op for `ms_confirmed_only=True`; `metric="ms_peptide_count"` uses
   regex `str.contains(allele)` against HLA strings with `*`, so literal
   alleles such as `HLA-A*02:01` do not match; and `metric="best_percentile"`
   silently falls back to allele-agnostic peptide counts if scoring imports
   fail. Filed as https://github.com/pirl-unc/tsarina/issues/34.

3. **Medium — `tsarina hits --allele` does not normalize user allele input.**
   `tsarina/cli_hits.py:271-275` post-filters by exact string. Current
   hitlist normalizes loaded `mhc_restriction` values, but tsarina does not
   normalize requested alleles, so `--allele A*02:01` drops canonical
   `HLA-A*02:01` rows. Filed as
   https://github.com/pirl-unc/tsarina/issues/33.

4. **Medium — legacy overlap helpers still treat scanner output as MS-only.**
   `tsarina/viral.py:416-442` and `tsarina/mutations.py:427-449` call
   `scan_public_ms()` and aggregate all hits without dropping
   `is_binding_assay`. This is the same class of issue as #31 but on older
   per-category helper APIs.

5. **Medium — local `hitlist` change has per-sample edge-case bugs.**
   `hitlist/pmhc_query.py:296-354` adds `query_by_samples()`, but a sample
   with an empty allele list is passed through to `query()`, where empty
   alleles mean "all alleles." The comment that empty sub-results remain
   visible to `groupby("sample_name")` is also false: inserting a scalar
   column into an empty frame still leaves zero rows. Filed as
   https://github.com/pirl-unc/hitlist/issues/188.

6. **Low — docs still carry stale pre-hitlist delegation details.**
   `README.md:121-125`, `README.md:292`, `docs/index.md:102-106`, and
   `docs/index.md:273` still mention `~/.tsarina` / `PERSEUS_DATA_DIR`
   and foreground `cancer_specific_viral_peptides()` where the current
   clinical default is stricter human-exclusive viral filtering. Filed as
   https://github.com/pirl-unc/tsarina/issues/32.

### Immediate Recommendation

Fix #31 first. It is foundational: `target_peptides()` feeds
`build_panel_matrix()`, the README's target-table example, and older target
selection APIs. The most direct fix is to mirror `personalize()` /
`export.build_ms_support_maps()`:

- default path: `indexing.load_ms_evidence(peptides=..., mhc_class=...,
  mhc_species="Homo sapiens")`
- explicit raw path: `scan_public_ms(...)`, then drop `is_binding_assay`
- keep source-flag aggregation shape stable for existing callers

Then fix #34 with targeted panel tests.

### Verification

- [x] `./format.sh`
- [x] `./lint.sh`
- [x] `./test.sh` — 214 passed

# Audit — Hitlist Alignment And TSA Goal

## Goal

Audit the current tsarina pipeline against two questions:

1. Do the recent hitlist-driven changes preserve the intended semantics of
   "MS-supported" evidence?
2. Does tsarina, as currently implemented, actually produce a
   high-confidence set of tumor-specific antigens suitable for vaccine and TCR
   therapy prioritization?

The audit should distinguish between:

- correctness of the code relative to the current hitlist/data interfaces
- scientific/product-fit of the resulting target set
- missing guardrails that could let non-tumor-specific or weakly supported
  candidates survive ranking/export

## Audit Plan

- [x] Review recent tsarina changes that touched hitlist integration, cached
      observations loading, public-MS scanning, and target ranking/export.
- [x] Trace the end-to-end pipeline from peptide generation through evidence
      collection, tissue-specificity filters, scoring, and target assembly.
- [x] Compare implemented filters and scoring rules against the stated goal of
      high-confidence MS-supported tumor-specific antigens for vaccines/TCRs.
- [x] Identify concrete failure modes, missing controls, and places where the
      current implementation could overstate confidence or tumor specificity.
- [x] Record findings with file/line references and recommended follow-up work.

## Verification

- [x] Read-only audit — no repo edits, so `format.sh`/`lint.sh`/`test.sh` N/A.

## Review

### Answers to the two audit questions

**Q1: Do recent hitlist-driven changes preserve the intended semantics of
"MS-supported" evidence?**  Yes. The refactor across v0.5.2 → v0.7.0 is clean:

- `is_binding_assay` filtering is consistent across fast path
  (`indexing.load_ms_evidence` drops binding rows by default) and slow path
  (`export.py:78-79`, `evidence.py:208-209`, `cli_hits.py:359-360` all filter
  after `scan()`).
- `--include-binding-assays` correctly routes to
  `hitlist.observations.load_all_evidence` on fast path (`cli_hits.py:367-376`)
  and skips the post-scan drop on slow path.
- `mhc_species="Homo sapiens"` is hardcoded in every library entry point
  (`personalize.py:232`, `export.py:86`, `evidence.py:205,216`); CLI exposes
  `--species`. No silent widening.
- `_species_kwargs`, `hla_only`, `human_only` shims fully removed (grep-clean).
- `resolve_cedar_path` silent-None on unregistered is by design.

**Q2: Does tsarina, as currently implemented, produce a high-confidence set of
tumor-specific antigens for vaccine / TCR prioritization?**  No — not via the
`personalize()` entry point.  The `targets()` path is closer but still has the
viral-exclusivity hole.  Several guardrails are missing, miswired, or defined
but unenforced.

---

### Findings — severity ordered, with file:line refs

#### CRITICAL (can let non-tumor-specific peptides into final output)

1. **`personalize.py:120,129` uses `cta_peptides()` instead of
   `cta_exclusive_peptides()`.** A 9-mer from MAGEA4 that also appears in a
   ubiquitously-expressed housekeeping protein is ranked as a CTA target.
   `targets.py:99,101` already uses the exclusive variant — the two paths
   diverge silently. `peptides.py:122-177` defines exclusivity as "peptide
   does not appear as a substring in any non-CTA protein-coding transcript".

2. **`personalize.py:177,180` AND `targets.py:117,121` both use
   `viral_peptides()` instead of `human_exclusive_viral_peptides()`.**
   `viral.py:245-302` defines the exclusive variant (filters viral k-mers
   that also match the human proteome); it exists and is documented but is
   never called from the main pipelines. An HPV16 E6 k-mer that collides
   with a human self-peptide can be ranked as a viral target.

3. **`personalize.py:266-288` — MHC scoring is effectively optional.**
   Scores are merged `how="left"` (line 286); `except ImportError: pass`
   (line 287-288) silently swallows an mhctools/mhcflurry install failure.
   A peptide with missing `presentation_percentile` then falls through to
   `line 303`: `(100 - percentile.fillna(100)).clip(lower=0) * 0.5`, which
   evaluates to 0 rather than excluding the row. The candidate still carries
   any MS-driven score and lands in the output with `best_allele=NaN`.

4. **`personalize.py:292-305` `priority_score` conflates modalities with
   uncalibrated weights and no tiering.**  A single MS hit anywhere earns
   `+10` (line 294); a peptide observed in cancer earns `+20` (line 296);
   a 0.5%-percentile prediction earns `~+49.75` (line 303-305). One MS hit
   plus `ms_in_cancer=True` (+30) beats a strong prediction with no public
   observation. There is no tier column, no minimum-confidence gate, and no
   distinction between "seen in one cell line" and "seen across many tissues".
   `scoring.py` defines `PRESENTATION_PERCENTILE_THRESHOLDS = (0.25, 0.50,
   1.0)` but nothing downstream consumes them.

#### MAJOR

5. **Healthy-tissue evidence is a soft penalty, not a gate
   (`personalize.py:297-298`: `-= ms_in_healthy_tissue * 50`).**  A peptide
   observed in normal tissue can still rank positively if MS-in-cancer,
   hit-count, and presentation contributions outweigh it. `targets.py:231-234`
   exposes a `cancer_specific` keyword that *does* hard-filter, but
   `personalize()` never wires one through.

6. **`mtec.py` and `negatives.healthy_tissue_peptides` are not wired into
   `personalize()`.**  mTEC is used in `selection.py:89,121` (sort-tiebreak)
   and joined by `export.py:127-133`, but the clinical `personalize()` →
   CSV path skips it. Thymic-self peptides (likely to be tolerized, so
   clinically useless as vaccine / TCR targets) are not flagged.
   `negatives.healthy_tissue_peptides` is defined but has zero callers
   outside its own docstring example.

7. **`personalize.py:307` — single-key sort, no deterministic tiebreak.**
   `sort_values("priority_score", ascending=False)`. Many candidates tie at
   `ms_hit_count * 10` values; pandas preserves insertion order, which
   depends on frame concat order in lines 115-194. Re-running the same patient
   inputs can reshuffle the top of the list.

8. **No per-CTA restriction / confidence gate in `personalize()`.**
   `personalize.py:122-126` admits any gene in `CTA_gene_names()` at
   `tpm >= min_cta_tpm`. Per-gene `restriction_level` / `restriction_confidence`
   computed in `tiers.py` are reporting-only in this path.  A LEAKY / LOW-
   confidence CTA contributes peptides indistinguishably from a STRICT / HIGH-
   confidence one.

9. **`min_cta_tpm=1.0` default (`personalize.py:50`) is below the
   `never_expressed < 2 nTPM` cutoff used to build the CTA CSV.**  Aligning
   runtime to 2.0 would close a small gap where a borderline-expressed CTA
   passes both filters.

#### MINOR

10. **Predictor-aware calibration is absent.**  `--predictor
    netmhcpan|netmhcpan_el` swaps the scorer, but `personalize.py:303` still
    treats percentile on a 0-100 scale with a fixed `*0.5` weight. NetMHCpan's
    percentile distribution differs from MHCflurry's; no warning is emitted.

11. **CEDAR silent-None** (`datasources.py:70-81`). Intentional per docstring,
    but a one-line stderr note when CEDAR is unregistered would help users
    distinguish "no CEDAR data" from "I forgot to register CEDAR".

#### POSITIVE (no action)

- `mutations.py:355-356` guards `mut_pep == wt_pep` → "mutant" peptides that
  are actually self are dropped.
- Fast-path ↔ slow-path parity for MS filters is clean.
- `_species_kwargs` / `hla_only` / `human_only` traces fully removed.
- `datasources.resolve_*` is the single registry entry point.

---

### Recommended follow-ups (grouped for separate PRs)

- **PR A (CRITICAL, 2 lines):** `personalize.py:120,129` → use
  `cta_exclusive_peptides`.  Add regression test that a peptide present in
  both MAGEA4 and a non-CTA protein does not appear in personalize output.
- **PR B (CRITICAL, ~6 lines):** `personalize.py:177` and `targets.py:117`
  → default to `human_exclusive_viral_peptides`. Add a kwarg
  `require_human_exclusive: bool = True` for explicitness.
- **PR C (CRITICAL):** Mandatory MHC-prediction gate in `personalize()`.
  Either drop NaN-percentile rows or annotate them with `tier="UNSCORED"`
  and exclude from top tiers. Fail loudly (not silently) on mhctools
  ImportError when `score_presentation=True` was requested.
- **PR D (CRITICAL + MAJOR):** Replace `priority_score` arithmetic with an
  explicit tier system consuming `scoring.PRESENTATION_PERCENTILE_THRESHOLDS`.
  Tier 1 = (MS cancer-only AND percentile ≤ 0.5) OR (≥3 MS hits AND ≤1%);
  Tier 2 = MS OR ≤2%; Tier 3 = else. Export a `tier` column. Deterministic
  sort `[tier, -ms_hit_count, presentation_percentile, peptide, best_allele]`.
- **PR E (MAJOR):** Add `enforce_tumor_specificity: bool = True` to
  `personalize()` that drops `ms_in_healthy_tissue=True` rows (mirrors
  `targets.py:231-234`'s `cancer_specific`).
- **PR F (MAJOR):** Wire `mtec.load_mtec_gene_table` and
  `negatives.healthy_tissue_peptides` into `personalize()` as flags-on-by-
  default filters; annotate excluded peptides in a diagnostic CSV.
- **PR G (MAJOR):** Add `restriction_level` / `restriction_confidence` gate
  for CTAs (default: require HIGH or MODERATE).
- **PR H (MINOR):** Raise `min_cta_tpm` default to 2.0; log CEDAR-missing
  once per invocation; emit a warning when `--predictor != mhcflurry`.

---

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

---

## Fix Panel Vital RNA Threshold (tsarina#43)

Goal: make `tsarina panel` use a biologically consistent default vital-tissue
RNA gate. The current `0.0` nTPM default treats HPA sub-1 nTPM background as a
hard veto, even though the CTA curation model uses deflated RNA and considers
somatic RNA detected at higher thresholds.

Plan:

- [x] Change the default `--vital-tissue-max-ntpm` / API default from `0.0`
      to `2.0`.
- [x] Update docstrings, CLI help, README/docs language, and tests so the
      documented default matches behavior.
- [x] Add focused tests showing sub-2 nTPM vital RNA does not exclude CTAs by
      default, while healthy-MS vital tissue evidence still gates non-allowlisted
      CTAs.
- [x] Bump version for the PR.
- [x] Run `./format.sh`, `./lint.sh`, and `./test.sh`.
- [ ] Open PR linked to tsarina#43; merge when CI is green and deploy.

Review:

- Local validation passed: `./format.sh`, `./lint.sh`, and `./test.sh`
  (`256 passed`). Live selector check confirms `NY-ESO-1` remains selected
  without the allowlist at the new 2.0 nTPM threshold, while MAGEA4 remains
  gated by gene-level healthy-MS heart evidence unless allowlisted.
