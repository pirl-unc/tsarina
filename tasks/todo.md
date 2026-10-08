# Plain-language documentation and linked navigation (2026-10-08)

Replace the homepage's "Start with your outcome" and "Continue to" copy with
literal workflow headings linked directly to their guides. Use short "See:"
references where a separate destination is useful. Review the README and
documentation introductions/headings for canned navigation, rhetorical
questions and vague narration; replace them with concrete descriptions.
Preserve command examples, scientific criteria, source tables and results.
Check any affected anchors and the rendered links before publishing.

- [x] Review the current documentation and create a feature branch.
- [x] Record the copy correction in lessons and file the documentation issue.
- [x] Edit navigation, introductions and headings; inspect the complete diff.
- [x] Run format/lint/full tests and strict MkDocs; verify rendered links/layout.
- [ ] Bump version, pass PR CI, merge, deploy and verify PyPI/Pages.

Plan check-in: linked workflow headings replace separate "Continue to" lines.
Plain descriptions replace canned lead-ins. Technical meaning stays intact.
Issue: https://github.com/pirl-unc/tsarina/issues/206.

Review: format/lint and strict MkDocs passed. The full real-model suite passed
684 tests with two optional skips and 19 warnings (85% coverage). All four
homepage workflow headings link directly to their guides; 134 internal links
and anchors across seven rendered pages and README relative destinations
resolve. The quoted navigation phrases are absent from docs and README.
Commands, tables, scientific criteria and result files are unchanged. Version
1.34.4 is prepared; final publication evidence will be recorded on the PR.

# Rendered HLA notation correction (2026-10-08)

The live evidence-comparison page parses paired asterisks inside HLA allele
names as Markdown emphasis, hiding required molecular notation. Use inline
code for every allele name on this page. Preserve all scientific results and
data hashes; verify the rendered text includes each exact allele string.
This is a documentation correction, with its own issue, versioned PR and
clean-main deployment.

- [x] Reproduce the live rendering error and create a feature branch.
- [x] File the bug, correct allele markup, and record the verification lesson.
- [x] Run format/lint/full tests and a strict MkDocs rendering check.
- [ ] Bump version, pass PR CI, merge, deploy and verify PyPI/Pages.

Plan check-in: the failure concerns rendered allele notation only. Inline
code preserves molecular names without changing analysis inputs or results.
Issue: https://github.com/pirl-unc/tsarina/issues/204.

Review: format/lint passed; the full real-model suite passed 684 tests with
two optional skips and 19 warnings (85% coverage). Strict MkDocs passed.
The rendered HTML preserves all eight molecular allele names exactly, and
all nine scientific assets retain their published source hashes. Version
1.34.3 is prepared; final merge/publication evidence will be posted to the PR.

# MS evidence, target HLA support and rare-CTA justification (2026-10-08)

Audit the published 2.5-kb designs without claiming that additional segments
fit the existing full construct. Quantify extra PRAME MS peptides, protein-level
HLA support and global allele support independently. For each retained protein,
report the change when it is removed from the final set, including mortality/
incidence expression bounds, unique MS peptides and the HLA carrier proxy.
Compare candidate segment allocations at explicit tolerated losses of modeled
expression coverage, preserving the current global allele set and RNA cap.
These allocation comparisons precede junction optimization and are not new
validated vaccine constructs. Preserve all biological/MS gates and source data.
Do not silently introduce a clinical benefit model or treat population carrier
probabilities as clinical protection. File the demonstrated allocator limits.

- [x] Inspect scoring precedence and audit additional PRAME evidence/HLA support.
- [x] Compute final-set removal contributions and top-three coverage baselines.
- [x] Compare coverage-loss tolerances and distinct-MS yield under the sequence cap.
- [x] Document per-protein justifications, assumptions and limitations; file issues.
- [x] Verify analysis independently, run required formatting/lint/tests, and review.
- [x] Publish the analysis through a versioned PR, PyPI and GitHub Pages.

Plan check-in: the current allocation counts distinct new MS peptides after
cancer gains and new global allele counts. Any positive cancer gain outranks
evidence regardless of magnitude. The current cancer score counts expression
once per protein and cannot represent extra HLA support specific to that target.
Compare alternatives before changing these scientific priorities.

Review: all eight allocation runs reached their model optimum. Preserve each
global allele's strongest restriction tier and current PRAME/MAGEA4/XAGE1A/B
target-specific support. At zero expression loss, strict keeps 93 observed
peptides with 19 proteins; loose reaches 89 with 18. A 0.05-pp tolerance gives
strict 11 proteins / 95 observed peptides / 49 PRAME peptides. These are native
allocations only, with no new assembled construct or junction validation.
Independent checks cover every candidate's source evidence, native coordinates,
coverage calculations, length accounting, 146/142 preserved restrictions per
variant, 98,500/98,312 non-CTA background occurrences, and 70,142 primary normal
8-mers. All pass. Two figures in PNG/SVG were visually checked; the docs page
and links were checked in Chrome. Format/lint, strict MkDocs and 684 full-suite
tests pass (two optional skips). Issues #201/#202 track the scientific allocation
limits; this analysis does not claim to resolve them. PR #203 merged and
version 1.34.2 was deployed from clean main; both published distributions
match the local release artifacts and PyPI SHA256 metadata. All seven PR
checks, main CI and the Pages workflow passed. Final live-page verification
found the allele-markup defect tracked in the correction above.

# Project aims and 2.5-kb vaccine budget (2026-10-07)

Use a literal project description in the shared website template and published
site. Keep the heading independent of design mode. The aim is one CTA vaccine
antigen prioritizing cancers responsible for the most deaths, with exact MS
ligand evidence and broad HLA support. Replace promotional slogans.

The full RNA cap is 2500 nt, including HBB 5′ UTR (50 nt), HBB_FI 3′ UTR
(268 nt), polyA (120 nt) and stop codon (3 nt). This leaves 686 encoded aa,
including the initiating methionine and linkers; maximum total is 2499 nt.
Regenerate all four comparisons with 686-aa/2500-nt caps using the audited
inputs and actual models, preserving source gates. Do not truncate an existing
construct. Refresh tables/plots and verify every changed artifact before PR,
clean-main PyPI release and Pages publication.

- [x] Replace headline/intro and capture the copy correction in lessons.
- [x] Refresh full-protein evidence with latest Hitlist 1.66.0.
- [x] Add compact-region tie-break without sacrificing higher coverage/evidence gains.
- [x] Regenerate strict/loose budget and ten-protein comparisons at the new caps.
- [x] Independently audit sequences, lengths, MS support and coverage tables.
- [x] Run format/lint/full tests, website checks and strict docs build.
- [x] Bump version, PR/CI/merge, deploy PyPI and verify live Pages.

Review: all four designs use 686 aa / 2499 total RNA nt with actual model
predictions and Hitlist 1.66.0 fresh full-protein evidence. Strict/loose budget
designs contain 23/26 proteoforms, 32/33 native pieces, 93/85 distinct MS peptides,
442/352 peptide–HLA pairs and 53/54 supported panel alleles respectively.
Equal-score allocation favors protein reuse and longer native pieces. PRAME
retains 139/57 aa in budget mode and 211/140 aa in ten-protein mode; the funnel
now separates piece omission from terminal trimming. Independent checks pass
for all four designs and all 175 coverage prefixes. Format/lint, 684 full-suite
tests (two optional skips), strict MkDocs, browser controls/layout and wheel/
sdist checks pass. Version 1.34.1 was merged and deployed; PyPI artifacts and all
361 published website files were verified. Release evidence:
https://github.com/pirl-unc/tsarina/pull/200#issuecomment-6048125061.

# Vaccine atlas and length-budget design (2026-10-07)

Specification: [vaccine-atlas-spec.md](vaccine-atlas-spec.md).

- [x] Regenerate strict/loose ten-target designs with current Hitlist; independently audit.
- [x] Verify plan: reusable report output plus standalone renderer and published example.
- [x] Implement coverage bounds, evidence/sample maps and length-budget allocation.
- [x] Package friendly website, source explanations, CT83 comparison and scientific plots.
- [x] Generate length-budget designs and independently reconcile results and figures.
- [x] Run format/lint/full tests, packaging, strict docs build and browser checks.
- [ ] Bump version, open PR, pass CI, merge and deploy from clean main to PyPI.
- [ ] Verify published GitHub Pages and PyPI artifacts; record review and next issue block.

Progress review: normal cardiac exclusions are verified against Atlas 2020.12
donor/sample rows and PMID 33858848. The default published examples use one
distinct donor and primary nonmalignant heart/brain/lung HLA-I evidence.
Both budget examples pass an independent 70,142-primary-8-mer check, including
their synthetic junctions. Strict retains 22 proteoforms / 94 exact MS peptides;
loose retains 32 / 133; both are 1000 aa / 3441 total RNA nt. CT83 passes
specificity and MS support under both definitions but loses the budget
allocation. All remaining predicted junction binders are explicitly reported.
The full model-enabled suite passes 683 tests with two optional skips.
Issues #193 (primary frequency values), #194 (longer normal ligands),
#195 (lazy figure layout), #196 (negative assay results) and #197 (legacy
binding-partition reads) are recorded and fixed in this branch. All four
designs, all 194 cumulative prefixes, cardiac exclusions and downloadable
artifacts reconcile independently. Strict MkDocs and browser-control checks
pass; an installed wheel renders a saved report without rerunning models.

# Restricted ten-proteoform vaccine design (2026-10-06)

Specification: [vaccine-supported-selection-spec.md](vaccine-supported-selection-spec.md).

- [x] Read current code, evidence audit, FDA reference and Hitlist issue status.
- [x] Create feature branch and specify exclusions, real-MS gate and target count.
- [x] Implement gene exclusions, explicit supported selection and whole-piece reservation.
- [x] Enforce positive MS modality and preserve accepted/rejected source audits.
- [x] Add regression tests and document the new selection contract/options.
- [x] Freshly scan expanded candidates and validate strict/loose ten-target designs.
- [x] Update figures/tables/examples and independently reconcile all scientific outputs.
- [ ] Run format/lint/full model tests, review, bump version and land PR.
- [ ] Deploy clean main, verify PyPI artifacts and record foundational follow-up work.

Plan check-in: keep MAGEA4 as the sole allowed MAGE-family target, preserve
CTA membership for excluded proteins, require genuine MS, and count ten
proteoforms that actually contribute native segments. Preserve the existing
1000-aa / 3500-nt limits; disclose feasibility failures rather than relax them.

Review: format/lint pass; 669 tests pass with actual MHCflurry, two optional
Topiary skips and 19 warnings (84% coverage). Hitlist 1.64.7's fresh scoped raw
scan covers 48,108 candidate peptides from the first 60 eligible positive-score
proteoforms per definition: 898 nonbinding records / 361 distinct peptides,
with positive MS modality subsequently enforced. Strict and loose each retain
ten distinct full-sequence targets in 20 native pieces, 1000 aa / 3441 nt.
Strict retains 113 distinct MS peptides / 639 pMHC assignments / 51 panel
alleles; loose retains 113 / 627 / 52. Counts are sequence/assignment units,
not independent patients. Strict includes CT83; loose instead includes CABYR.
Only MAGEA4 remains eligible among MAGE genes.

Independent checks pass for all 34 runtime artifact hashes per definition,
translation/UTR/stop/polyA lengths, target count, MS modality and sample
assignments, retained ligand coordinates, all final junction windows and
independent scans of 98,500 non-CTA translated occurrences. All six scientific
figures have been visually checked. Junction predictions below 1000 nM decline
1068→517 (strict) and 1232→656 (loose); both explicitly require review.
Bundled reports include cancer incidence/mortality/p95, complete sequence
funnels, all candidate outcomes, per-allele counts, evidence assignments,
accepted/rejected observations, figures and FASTAs. Hitlist #644 remains
upstream; this PR closes Tsarina #189 through an explicit consumer MS gate.
Python 3.9 syntax, strict MkDocs, distribution build/Twine checks and all
48 bundled snapshot hashes also pass. The sdist contains exact copies of all
report artifacts. Merge and clean-main deployment remain pending; final
release evidence will be recorded on the PR without a direct main commit.

---

# CTA vaccine design (2026-10-06)

Detailed specification: [cta-vaccine-spec.md](cta-vaccine-spec.md).

- [x] Inspect OncoRef, Pirlygenes, Tsarina, Vaxrank, scientific predictors and lessons.
- [x] Create feature branch and document objective, identity, funnel and constraints.
- [x] Implement mortality/proteoform input adapter and strict/loose definitions.
- [x] Implement native interval subtraction and auditable panel MS support.
- [x] Implement padding/order/linker optimization, cleavage and final audit.
- [x] Implement DNA/RNA assembly, constraints, CLI and reports/figures.
- [x] Add regression coverage, offline example and real scientific validation.
- [x] Run format/lint/full tests, review and bump to 1.33.0.
- [ ] Open PR, pass CI, merge, deploy and verify published artifacts.
- [ ] Record review evidence and dependency-ordered follow-up work.

Plan check-in: use OncoRef as the CTA and expression authority; retain native
CTA-exclusive intervals and measured MS provenance; audit actual assembled
sequence and explicit limitations. No gene-level percentile arithmetic.

Review: PR #188 adds the feature and bumps 1.32.2 to 1.33.0. The compatible
project environment passes pip check with imported Hitlist 1.64.7, OncoRef
1.8.207, PyEnsembl 2.23.2 and mhctools 3.47.1. Format/lint pass; the final
model-enabled suite passes 654 tests (two optional Topiary backend skips,
19 warnings; 84% coverage). Regression checks also prevent clamped padding
from filling the beam with duplicate constructs and allow linker/padding
changes together at the length cap. Python 3.9 syntax, synthetic artifacts
and figure QA pass.

Same-gene HSCHR annotations no longer falsely veto PRAME/MAGEA3/MAGEA6 in
this feature; expression keys and independent non-CTA vetoes remain intact
(legacy partition follow-up #187). A full current-index rebuild completed both
raw scans and mappings, but its lossless contributor export exceeded about
58 GB available scratch space (Hitlist #643). It was stopped before filling
the disk; only its owned temporary graph was removed. The stale observation
cache was neither read as current nor restamped.

Final real-data API validation instead uses explicit fresh scoped Hitlist
1.64.7 scanner/collector input: 588 human class-I MS rows, 217 peptides,
1,101 contributor links, source hashes, current classification/dedup/exclusion
contracts. Strict and loose both retain eight of ten ranked proteoforms in
26 native pieces: 1000 aa / 3441 nt with HBB/HBB_FI and 120-nt polyA.
Remaining below-1000 nM window/allele predictions are 1088 (from 1481) and
950 (from 1444), respectively; both designs are explicitly review-required.
GAGE1/GAGE2A lack qualifying MS support; no automatic backfill occurs.
Independent checks verified translation, native/assembled ligand coordinates,
all final window/allele pairs, constraints and 60 output hashes. A separate
scan of 98,500 independent translated background occurrences confirms no
retained native 8-mer collision. Bundled docs include complete selected-cancer
tables, funnels, layers, ligand/junction/cleavage tables, FASTAs, figures and
hashed validation snapshots; full provenance bundles remain under ignored
vaccine-designs/2026-10-06/. Merge/PyPI verification remains pending.

---

# Compatible development installs (#184, 2026-10-05)

## Specification

Resolve Tsarina with its dev extra and every present scientific-data sibling
as one editable installation. Before changing installed packages, use pip's
dry-run report to project the complete resulting environment, including
untouched installed consumers. Validate that projection with pip check in a
temporary metadata-only virtualenv, so a stale sibling cannot violate an
existing consumer's minimum unnoticed. Reject conflicts with the consumer,
requirement, proposed version, and instructions to update the checkout or use
a separate virtualenv. Install the checked editable set together with normal
dependency resolution; require a successful final pip check. Use the selected
Python interpreter consistently and pip >=23.0's stable report format.

Retain active-virtualenv selection, absent-sibling release fallback, optional
predictor policy, and resolved-import reporting. Avoid pip internals, custom
requirement parsing, runtime dependency additions, and individual --no-deps
installs. A failed preflight must preserve installed distributions. Regression
tests will exercise actual pip resolution offline in disposable virtualenvs,
including stale sibling-vs-sibling and installed-consumer constraints, newly
added dependencies, editable preservation, and final-check failure handling.
Release as 1.32.2 via PR and deploy from clean main.

## Plan

- [x] Read #184, local workflow/lessons, installer APIs, and current environment.
- [x] Implement joint resolution and whole-environment preflight before install.
- [x] Prove rejection and successful editable/dependency installs with regressions.
- [x] Run format, lint, full tests including real models, and review the diff.
- [ ] Bump version, open PR, wait for CI, merge, and deploy clean main.
- [ ] Verify PyPI artifacts/issue closure and review dependency-ordered next work.

## Plan check-in

The current environment passes pip check. Joint sibling resolution alone
cannot protect an untouched installed consumer, and installer success alone
does not establish a consistent environment. The implementation will use
pip's public report, interpreter-selection, and check commands for those
two separate gates, with no preflight package installation.

## Review

Initial offline regressions exposed a local interpreter constraint: copying
uv's macOS Python executable loses its relative libpython dylib. Use the
POSIX symlink behavior of `python -m venv` in disposable environments; the
resolution/preflight design is unchanged. Verification continues with that
correction before any real development installation.

The three core regressions fail against main's original script: stale sibling
installs exit successfully after downgrading mhcgnomes, and newly declared
dependencies remain absent. All six initial tests pass on the repaired script.
The first full format/lint/model-enabled test run passed 615 tests (19 warnings).
The real sibling preflight refused mhcgnomes 3.64.2 against hitlist 1.64.6's
>=3.64.4 minimum, preserved every installed distribution version, and left
pip check passing. No sibling was updated to disguise that expected refusal.

Review also isolates the metadata-check subprocess from PYTHONPATH/PYTHONHOME
so foreign distribution metadata cannot shadow proposed versions; an added
negative regression exercises that case. Version 1.32.2 is bumped explicitly
on the feature branch because the current deploy.sh uploads the existing
version (the script's behavior is already documented in this task log).
Final validation, PR, merge, and release verification remain pending.

Final validation: format/lint pass; the full model-enabled suite passes all
616 tests (19 warnings, 97.39 seconds; 81% runtime coverage). Seven new tests
cover the two rejection paths, successful installation of a new dependency
while preserving all editables, absent-sibling fallback and environment
markers, final-check failure, future report rejection, and PYTHONPATH metadata
shadowing. The helper is included in local/CI lint gates. Python 3.9 syntax
checks and Bash syntax checks pass. Reviewed the complete change: no runtime
dependency or scientific-selection behavior changes; normal installation
resolves the joint graph and constrains it to the preflight-checked versions.

Built wheel and sdist successfully before the PR. The sdist includes byte-for-
byte copies of develop.sh, scripts/develop.py, and its regression tests; wheel
version matches 1.32.2. Post-merge clean-main deployment, artifact digests,
issue closure, and next-work review will be recorded in the PR description,
so release evidence does not require a direct commit on main.

# Correct fetch-all destination output (#180, 2026-09-29)

## Specification

Keep hitlist as the owner of download paths and its existing progress output.
`tsarina data fetch-all` must not describe the independently configured corpus
directory as the asset-cache destination. Remove the directory from tsarina's
redundant summary, retaining the number of fetched files (including zero).
Preserve force forwarding and download error propagation. No cache migration,
directory creation, dependency change, or new download API is required.
Release the fix as 1.32.1 through a PR, then deploy from clean main.

## Plan

- [x] Inspect #180, current main, upstream download output, and release scripts.
- [x] Correct the summary and verify distinct corpus/cache paths without downloads.
- [x] Run format, lint, and full tests with real-model checks and one worker.
- [ ] Review the diff, merge after CI, deploy, and verify PyPI artifacts.
- [ ] Confirm issue closure and inspect the next relevant open work.

## Review

Verified the real upstream progress formatter through the tsarina parser with
a temporary asset and an independently configured corpus directory. Both force
modes preserve the actual asset-cache output without consulting `data_dir()`.
An empty asset registry prints a zero count, and download errors still propagate.
Format and lint pass. `TEST_SH_MAX=1 ./test.sh --run-mhcflurry` passed all 609
tests, including the real-model contracts (19 warnings, 70.42 seconds).
Reviewed the complete diff: the runtime change removes only the incorrect
directory suffix, and no new dependency API is needed. Completion claims
identify the released issue and version explicitly; PR #183 did not include
this fix. Merge, publication, and artifact verification are tracked in PR #185.

# Review fixes and release (2026-09-24)

## Specification

Resolve #175–#178 together in a versioned PR from a feature branch. Screen
mutant targets against all reference-human coding proteins before they can
be ranked as non-self; retain an explicitly documented raw-enumeration mode
for audit use. Reuse the existing candidate-driven overlap scan rather than
building a full multi-length proteome posting index. Verify PLK1 and RAS-family
overlaps against Ensembl 112 and show that patient/unified target paths use
the safe default.

Give `build_panel_matrix` an explicit, validated maximum presentation percentile
for binder counts, defaulting to the existing candidate threshold of 1.0.
Exclude missing/non-finite scores and test threshold boundaries, zero-count
alleles, duplicates, and preservation of the other metrics.

Replace the invalid NetMHCpanEL import with the supported upstream predictor
API while preserving presentation/affinity meanings. Exercise the advertised
selectors with actual installed libraries and, when available, the binary.
Trace the single failing MHCflurry calibration test: distinguish the public
2.2.1 API from the development checkout, retain real calibration verification,
and make the optional integration contract explicit rather than simply skipping
the assertion. File independently reproduced upstream library bugs with their
owning repos and link them in the PR.

## Plan

- [x] Reinspect current files, prior evidence, lessons, and dependency APIs.
- [x] Finalize dependency/test contracts and capture failing regression cases.
- [x] Implement human-overlap filtering, binder counts, and predictor fixes.
- [x] Repair the calibration test contract and document supported behavior.
- [x] Run format, lint, full tests, real-data checks, and review the full diff.
- [ ] Bump version, create/review PR, merge after CI, deploy from clean main.
- [ ] Verify PyPI artifacts and issue closure; identify next upstream work.

## Review

The reported baseline failure was
`tests/test_alleles.py::test_global53_default_uses_mhcflurry_runtime_calibration_when_available`,
an AttributeError for `Class1AffinityPredictor.percent_rank_calibrated_allele`
on installed MHCflurry 2.2.1. Tracked explicitly as #178. Direct calls to the
public `percentile_ranks()` API verified every default allele and calibration
reuse on that same installed release. The test now exercises those numeric
results, checks equivalent C*14 ranks, and still rejects uncalibrated C*15:05.

Fourteen new mutation/filter regression cases failed against the original
implementation. The shared coding-protein scan now removes exactly the 14
matches found by the independent full-FASTA audit (722 raw hotspot rows ->
708 screened rows). Tests cover alternate isoforms, IG segments, explicit raw
enumeration, reference errors, and the patient/unified pipeline boundaries.
An initial test-fixture context omitted the residue before the BRAF window;
corrected it to the actual Ensembl sequence before final verification.

The NetMHCpan selector defect belongs to Tsarina: no generic `NetMHCpanEL`
export was found in the inspected upstream history. Corrected #177's title
instead of asserting an unverified upstream removal. Both selectors now use
the working version-detecting adapter, with EL and BA metrics kept distinct.
The MHCflurry failure also belongs to the downstream test, not to the library.

`pip check` reports unrelated pre-existing shared-environment conflicts
(eureka-bench pins, dataclasses-json/marshmallow, TensorFlow/h5py, tweety,
pdfx). These are installed-version conflicts, not demonstrated upstream code
failures; no dependencies were replaced to mask them. Previously tracked
sercol#4 and topiary#374 are now closed and those conflicts no longer appear.

Final local verification: `./format.sh` and `./lint.sh` pass; full suite is
**605 passed, 19 warnings, no skips**. Both real NetMHCpan-selector integration
cases and the previously failing MHCflurry calibration test pass. A separate
real MHCflurry 2.2.1 smoke prediction returns numeric scores for SLYNTVATL and
poly-A. Wheel/sdist build and `twine check` pass for 1.31.7. Full-diff review
confirmed that existing CTA scan policy and other matrix metrics are preserved.

# Major-issue review (2026-09-24)

## Specification

Review current main (e71c1b2, 1.31.6) for consequential correctness,
integration, and release problems. Focus on candidate selection, evidence
assignment, safety filters, and deployment/install contracts. Validate each
finding with a minimal reproduction or an unambiguous execution path;
distinguish current defects from historical issues and environment drift.
This is an audit, with no requested product-code changes or release.

## Plan

- [x] Inspect repository state, instructions, lessons, and recent changes.
- [x] Establish lint/test baseline and inspect resolved dependency paths.
- [x] Trace high-impact selection/evidence/scoring paths and boundary cases.
- [x] Reproduce major findings; check existing issues and file new defects.
- [x] Record review evidence, limitations, and prioritized findings.

## Review

Reviewed current code rather than carrying forward historical findings. No
product-code edits, PR, or release are part of this audit. Findings are open:

1. **P1 — mutant reference-human overlap**
   ([#175](https://github.com/pirl-unc/tsarina/issues/175)). Mutant generation
   compares only against the source transcript's wild-type k-mer. Searching
   all 722 generated hotspot rows against 123,495 Ensembl 112 reference
   proteins found 14 overlapping rows: KRAS G12R in RHOT2/RASL10B and BRAF
   V600K in PLK1. With public MS simulated as absent, actual peptide generation
   and NetMHCpan scoring retained the PLK1-identical `KIGDFGLATK` as STRONG
   for HLA-A*03:01 (presentation percentile 0.027). This demonstrates sequence
   non-exclusivity, not healthy-tissue presentation or clinical toxicity.
2. **P2 — binder-count matrix has no binding gate**
   ([#176](https://github.com/pirl-unc/tsarina/issues/176)). The older public
   `build_panel_matrix(metric='peptide_count')` counts all predictions,
   including weak/unscored pairs. Reproduced equal counts of two for strong
   A*02:01 and 95th/99th-percentile B*07:02 predictions. This is distinct from
   the main `tsarina panel` command and from closed issue #34.
3. **P2 — advertised EL predictor adapter is incompatible with mhctools**
   ([#177](https://github.com/pirl-unc/tsarina/issues/177)). Installed mhctools
   3.44.55 lacks `NetMHCpanEL`. Real `netmhcpan_el` scoring fails at import;
   the same input succeeds with `netmhcpan`, including presentation and
   affinity output, proving the executable is installed and working.
4. **Baseline test/dependency contract failure**
   ([#178](https://github.com/pirl-unc/tsarina/issues/178)). Full
   `TEST_SH_MAX=2 ./test.sh`: 585 passed, one failed, 19 warnings. The optional
   calibration test expects `percent_rank_calibrated_allele`, absent from
   installed MHCflurry 2.2.1. The development checkout has the method; this
   is not evidence of a target-selection failure. No environment packages
   were replaced to hide the baseline failure.

Lint passed. All 19 hotspot transcript reference residues also matched the
cached Ensembl 112 proteins. Dependency paths were inspected: tsarina and
hitlist resolve locally; oncoref 1.8.194, mhcgnomes 3.64.1, pyensembl 2.10.4,
and MHCflurry 2.2.1 resolve to installed packages. Reproduction scripts used
real sequence data and real NetMHCpan where stated; panel-count scores were
controlled fixtures. No claim is made that this was an exhaustive audit or
that the entire suite passed. Fixes and release verification remain future
work tracked in the linked issues.

# Issue release series (2026-09-22)

## Specification

Create, review, merge, and publish one PR for each open issue: #132, #146,
#160, #130, #131, and #120. Audit current behavior before implementing an old
proposal: oncoref now owns CTA definitions/proteoforms, and #161 already added
flagged clinical targets. Each PR must contain a version bump and pass
`./format.sh`, `./lint.sh`, and `./test.sh`; deploy with `./deploy.sh` from
clean main after merge and verify both PyPI distributions.

Use current development checkouts for every locally developed dependency,
including optional predictors when installed. Compare local branches with
their remote heads, preserve unrelated work, and record resolved imports.
File newly discovered defects in their owning repositories and link them in
the relevant PR. Read primary literature for biological interpretation;
preserve evidence provenance and distinguish candidates from approved sets.

## Plan

- [x] Audit and update local dependency environment; establish full-suite baseline.
- [x] #132: verify oncoref ownership, dependency minimum, required integration
      coverage, and migration guidance; PR, review, merge, deploy.
- [x] #146: update Actions to supported runtimes across workflows; verify
      the Python 3.9–3.12 CI matrix; PR, review, merge, deploy.
- [x] #160: audit all CTA drop paths and flagged-target caveat propagation;
      reproduce remaining gaps, fix and test; PR, review, merge, deploy.
- [x] #130: audit affected genes against current oncoref protein evidence;
      correct stale upstream inputs if present, verify downstream parity;
      PR, review, merge, deploy.
- [x] #131: check SUN5/SUN3/SPAG4 evidence and literature, correct candidate
      coverage at its owner while retaining somatic caveats; PR, review,
      merge, deploy.
- [x] #120: document HERV-K locus/family coverage limits and appropriate
      evidence sources from primary literature; open and review PR #174.
- [ ] Merge #174, deploy 1.31.6, and record the post-merge publication audit
      in [PR #174](https://github.com/pirl-unc/tsarina/pull/174).
- [ ] Audit all six issue/PR/release states and identify next dependent work.

## Review log

### #120 implementation and verification plan

- [x] Inspect CTA and viral generation: no locus/family HERV-K quantifier or
      adapter exists; single gene identifiers do not cover the HML-2 family.
- [x] Read Telescope and ERVmap methods, and HERV-K experimental antigen
      evidence. Keep RNA abundance, translated ORF, and presentation distinct.
- [x] Document the explicit current coverage limit at curation and workflow
      entry points, the appropriate separate evidence-source boundary, and
      the provenance needed by a future integration. Do not imply an adapter
      is implemented or that gene-panel absence is negative ERV evidence.
- [x] Bump to 1.31.6; build the docs, run format/lint/full tests, review the
      claims against sources, and open separate PR #174.
- [ ] Merge/deploy and attach the final release verification to PR #174.

### #131 implementation and verification plan

- [x] Audit current SUN-domain coverage: SUN5 is strict TESTIS, SUN3 is
      retained but excluded with somatic IHC; neither requires another seed.
- [x] Review human SPAG4 evidence (PMIDs 14614621, 23602831), SUN5 colorectal
      evidence (PMID 36358787), and HPA v23. Human SPAG4 pancreas expression
      rules out strict testis selection; do not extrapolate mouse specificity.
- [x] File oncoref#548 for SPAG4 candidate-reference discoverability, retaining
      its normal-pancreas caveat without promoting it into a default set.
- [x] Add regression coverage for SUN5/SUN3 presence, canonical membership,
      Ensembl identity, and distinct RNA/protein evidence; document decisions.
- [x] Bump to 1.31.5, run format/lint/full tests, review, open the separate PR,
      then merge/deploy and verify artifacts.

### #130 implementation and verification plan

- [x] Recompute the two historical stale protein columns for all ten genes
      from oncoref's pinned HPA v23 normal-tissue IHC. Every current value
      matches: no upstream data repair is required.
- [x] Add downstream regression coverage for real IHC availability and
      canonical evidence propagation on each affected gene. Keep membership
      decisions separate from protein evidence; do not promote excluded genes.
- [x] Document the current owner, audit result, and protein/membership
      distinction. Verify the supported oncoref floor carries the corrected data.
- [x] Bump 1.31.4, run format/lint/full tests, review the diff, open a separate
      PR, merge/deploy and verify publication.

### #160 implementation and verification plan

- [x] Reproduce #169: a strict CTA present in the clinical registry returns
      through the flagged path after failing confidence or mTEC selection.
- [x] Classify strict versus clinical-only genes before selection. Preserve
      explicit clinical-only requests with visible oncoref/overlap caveats;
      apply the optional mTEC gate and tumor TPM gate to both categories.
- [x] Emit specific warnings for unknown symbols, upstream exclusions,
      absent default expression, restriction confidence, mTEC, and low TPM.
- [x] Verify mixed inputs, missing TPM, genuine flagged targets, and no
      peptide generation for dropped genes. Correct the documentation's
      afami-cel target claim against primary clinical sources.
- [x] Bump to 1.31.3; run format/lint/full tests; review and open a separate
      PR closing #160 and #169; merge/deploy after release authorization.

### #146 implementation and verification plan

- [x] Upgrade checkout/setup-python in all three workflows to their current
      Node.js 24 releases (v7); upgrade Pages upload/deploy actions to v5 so
      documentation publishing also avoids deprecated nested runtimes.
- [x] Preserve triggers, permissions, commands, and Python 3.9–3.12 coverage.
- [x] Validate YAML and each action's published runtime/compatibility contract.
- [x] Bump 1.31.1 → 1.31.2; run format, lint, and the full development suite.
- [x] Review the complete diff, open a separate PR, verify all CI jobs, merge,
      deploy from clean main, and verify both PyPI artifacts.

- Initial state: clean main at 3c3a5d8, version 1.31.0; all six issues open.
- GitHub CLI works with network/keychain access outside the sandbox.
- Current `deploy.sh` builds/uploads the existing version; it does not bump
  or commit despite AGENTS.md's older description. Bumps belong in each PR.
- #132 review: migration code already landed in #142/#143. Added explicit
  migration guidance, obsolete-sync replacement, and executable integration
  checks; bumped 1.31.0 → 1.31.1. No runtime duplication is needed.
- Verified remote/current branches and imports for hitlist, oncoref, mhcgnomes,
  pyensembl, datacache, gtfparse, serializable, sercol, mhcflurry, mhctools,
  topiary, varcode, and osteosarc. All resolve to local development checkouts.
- Dependency metadata problems are tracked in openvax/sercol#4 (existing
  serializable/simplejson constraints) and openvax/topiary#374 (new report:
  exact osteosarc pin replaces its current development checkout). Editable
  sources remain installed; these unresolved metadata constraints are not
  represented as a clean dependency-resolution check.
- #132 gates passed: format, lint, full suite (561 passed, 7 warnings),
  including the optional live mhcflurry calibration test. Reviewed the diff
  against #132 and the current ownership implementation; no remaining code gap.

---

# PR — Reorganize Documentation from Overview to Reference (2026-07-24)

## Goal

Make the documentation readable in progressive layers: first explain what
Tsarina does and which workflow a reader should choose, then present operational
guidance, and only then expose schemas, flags, formulas, and edge cases.

## Audit findings

- Tracking issue:
  [tsarina #144](https://github.com/pirl-unc/tsarina/issues/144).
- `README.md` and `docs/index.md` duplicate nearly the same 400-line manual and
  have already drifted in product naming and CLI examples.
- The largest block, CTA × HLA panel internals, is embedded in both entry
  documents instead of a focused workflow guide.
- Workflow choice is implicit; users encounter target taxonomy before learning
  whether they need personalized target selection or cohort panel design.
- Data/evidence concepts, scoring, naming, and development reference are
  separate top-level fragments rather than one coherent reference layer.
- `docs/curation.md` is accurate but moves directly into implementation
  ownership without an at-a-glance boundary and reader-oriented purpose.

## Plan

- [x] Make `README.md` a concise product overview with a workflow map, quick
      start, target-category summary, documentation map, and development entry.
- [x] Make `docs/index.md` a reader-oriented documentation hub: choose a
      workflow first, understand the shared pipeline second, then navigate to
      focused guides.
- [x] Create a personalized-target guide organized as outcome → required
      inputs → basic workflow → output/prioritization → advanced considerations.
- [x] Create a panel-design guide organized as outcome → default pipeline →
      selection stages → output controls → evidence tiers → HLA/coverage
      reference → advanced options.
- [x] Create a data-and-evidence reference organized as evidence model → source
      classification → data setup → tissue/scoring helpers → output naming.
- [x] Reorganize CTA ownership guidance as summary → ownership boundary →
      evidence flow → API behavior → maintenance.
- [x] Verify documented commands and flags against current CLI help, validate
      Markdown links/headings, and remove stale duplicated wording.
- [x] Bump the package version and run `./format.sh`, `./lint.sh`, and
      `./test.sh`.
- [ ] Open, merge, and deploy the PR.

## Review

- Replaced the duplicated 400-line README/index manuals with a 107-line
  product overview and a 120-line workflow-oriented documentation hub.
- Moved operational detail into focused personalized-target, panel-design, and
  data/evidence guides, each ordered from purpose and defaults to advanced
  reference material.
- Reordered CTA curation guidance around an at-a-glance ownership boundary
  before its schema, API, and maintenance details.
- Corrected stale product naming, personalized output columns, panel progress
  flags, the data discovery command, inclusive panel cutoffs, and the hotspot
  registry's gene count.
- Verified all local Markdown links and every documented CLI option against the
  current parsers.
- Released code version is prepared as 1.24.1.
- Verification passed: `./format.sh`, `./lint.sh`, and `./test.sh` (423 tests).

---

# PR — Make oncoref the Sole CTA Definition Authority (2026-07-24)

## Goal

Tsarina must consume cancer-testis antigen definitions from current `oncoref`
instead of maintaining a second CTA curation implementation or bundled copy.
Tsarina may retain only target-selection evidence that belongs to its own layer,
such as healthy-tissue immunopeptidomics annotations keyed by Ensembl gene ID.

## Design

- Raise the runtime floor to the current released `oncoref` API and import CTA
  and proteoform functions from their semantic submodules.
- Delegate canonical/default/filtered/unfiltered/excluded/low-expression CTA
  membership and alias resolution directly to `oncoref.cta`.
- Keep Tsarina-specific axis queries and MS-aware evidence enrichment, but
  restrict canonical membership with IDs returned by `oncoref` rather than
  reconstructing oncoref's specificity decisions.
- Replace the bundled full CTA/HPA curation table with a narrow MS-evidence
  overlay containing only `Ensembl_Gene_ID` and `ms_*` columns.
- Remove obsolete local CTA curation/regeneration scripts and constants whose
  authority has moved to `oncoref`.
- Remove the Tsarina-only H1-6 CTA evidence row. If its exclusion evidence is
  absent upstream, file an `oncoref` issue instead of preserving a competing
  local CTA universe.
- Update tests and documentation so ownership, counts, and downstream
  enrichment boundaries are explicit.

## Verification

- [x] Search open `oncoref` issues for any discovered upstream gaps; filed
      [oncoref #435](https://github.com/pirl-unc/oncoref/issues/435) for the
      reproducible restriction-confidence synthesis error.
- [x] Assert every foundational Tsarina CTA helper exactly delegates to
      `oncoref` and that Tsarina evidence has the same CTA row universe.
- [x] Assert the MS overlay has a narrow schema and unique unversioned gene IDs.
- [x] Assert the packaged wheel contains no duplicate full CTA curation table.
- [x] Run focused integration tests against oncoref 1.8.150 (131 passed).
- [x] Run `./format.sh`.
- [x] Run `./lint.sh`.
- [x] Run `./test.sh` (423 passed).

## Review

- `oncoref>=1.8.150` now owns the CTA row universe, all membership tiers,
  aliases, HPA axes, tissue definitions, and proteoform groups.
- Tsarina ships only a four-column generic MS overlay and HPA cancer-prevalence
  feature tables with no membership or specificity fields.
- The wheel contains no CTA-definition or proteoform fallback table.
- H1-6 retains generic MS evidence but cannot re-enter the CTA universe.
- The audit found no missing upstream CTA row: the seven local-only rows were
  intentional oncoref histone/tubulin-family exclusions.
- oncoref #435 tracks the discovered error where reproductive RNA could
  incorrectly boost a SOMATIC protein restriction call.

---

# PR - Adopt Oncoref CTA Defaults (2026-07-10)

## Goal

Make tsarina's default CTA gene set agree with `oncoref.cta_gene_ids()` while
retaining tsarina's mass-spec safety evidence. CSH1 and H1-6 must remain
inspectable as non-default candidates with their supporting safety evidence, but
neither should be returned by the default CTA panel helpers.

## Spec

- Treat oncoref's CTA evidence and specificity decisions as canonical for
  default CTA inclusion.
- Preserve tsarina's `ms_*` evidence columns by joining the local bundled
  evidence onto the oncoref CTA frame.
- Preserve H1-6 as a tsarina-only excluded evidence row so its
  `RECURRENT_HEALTHY` mass-spec signal is not lost even though oncoref excludes
  histones from the canonical CTA universe.
- Ensure CSH1 follows oncoref's exclusion (`passes_filters=False`,
  `specificity_status=excluded_normal_expression`) while keeping its lung and
  smooth-muscle RNA safety evidence visible.
- Update public accessors, automatic panel filtering, docs, and tests so the
  default set is 293 genes and matches oncoref exactly.

## Plan

- [x] Add oncoref as a runtime dependency and load oncoref CTA evidence as the
      base CTA frame.
- [x] Join tsarina-local `ms_restriction`, `ms_healthy_somatic_tissues`, and
      `ms_pmids` columns onto the oncoref frame.
- [x] Append H1-6 from tsarina's local bundled evidence as an excluded
      tsarina-only candidate with specificity metadata.
- [x] Add canonical CTA mask helpers and use them for default/filtered CTA
      accessors and automatic panel candidate filtering.
- [x] Update tests for oncoref count parity, CSH1/H1-6 exclusion, and retained
      mass-spec evidence.
- [x] Update docs/readme wording for oncoref-backed CTA curation and counts.
- [x] Rename the proteoform registry sync/test authority to oncoref so the
      mirror parity test runs in this PR.
- [x] Run `./format.sh`, `./lint.sh`, and `./test.sh`.

## Review

- `CTA_gene_names()` and `CTA_gene_ids()` now exactly match oncoref's canonical
  default sets (293 genes).
- `CTA_filtered_gene_names()` / IDs now match oncoref's canonical filtered tier
  (302 genes), while `CTA_evidence()` keeps a 391-row evidence universe with
  tsarina's MS safety columns.
- CSH1 is excluded through oncoref's `excluded_normal_expression` decision and
  keeps its lung/smooth-muscle RNA safety evidence visible.
- H1-6 is appended as a tsarina-only excluded evidence candidate with
  `ms_restriction=RECURRENT_HEALTHY`, 15 healthy somatic tissues, and 48 PMIDs.
- Automatic panel selection and CTA partitioning use the canonical default mask
  rather than raw `passes_filters`.
- Verification passed: `./format.sh`, `./lint.sh`, and `./test.sh`
  (474 passed). The proteoform registry parity test now uses oncoref, so it
  runs under the PR dependency set instead of being skipped.

---

# PR - Deconvolve Multi-Allelic MS Evidence (2026-05-07)

## Goal

Make all multi-allelic MS evidence prediction-deconvolved rather than treating
every listed allele as equally supported. The ranking should be
monoallelic/exact MS > multi-allelic sample/restriction deconvolution >
unrestricted class-I MS prediction > pure prediction.

## Plan

- [x] Score exact multi-allele restriction rows against every listed allele
      needed to choose a best-of-restriction allele.
- [x] Attribute non-monoallelic exact restriction rows only to the best
      predicted listed allele, not every listed panel allele.
- [x] Preserve current sample-narrowed donor-bag behavior and unrestricted
      class-I fallback.
- [x] Add regression coverage for exact multi-allele rows where multiple panel
      alleles pass the sample cutoff.
- [x] Update docs/CLI wording to describe the evidence ladder.
- [x] Bump the patch version.
- [x] Run ``./format.sh``, ``./lint.sh``, and ``./test.sh``.
- [ ] Open, merge, and deploy the PR.

## Review

- Exact multi-allele restriction rows now score all listed HLA alleles needed
  for deconvolution when the row overlaps the selected panel.
- Non-monoallelic exact rows now assign ``sample_allele_ms`` only to the best
  predicted listed allele rather than every listed panel allele.
- Sample-narrowed donor-bag rows keep the same best-of-set behavior, and
  class-only rows with no usable exact/donor set stay ``unrestricted_ms``.
- Docs and CLI help now describe the intended tier ladder:
  monoallelic MS > multi-allelic deconvolved MS > unrestricted class-I MS >
  pure prediction.
- Verification passed: ``./format.sh``, ``./lint.sh``, and ``./test.sh``
  (307 tests).

---

# PR - Handle Peptide Attribution Provenance (2026-05-07)

## Goal

Recover sample-narrowed hitlist observations whose provenance is
``peptide_attribution`` and make the sample-bag deconvolution rule explicit:
donor allele bags are candidate sets, and selected pMHC evidence comes from
prediction against the bag, not from treating the bag as an HLA assignment.

## Plan

- [x] Treat ``peptide_attribution`` rows like ``sample_allele_match`` rows when
      extending the allele scoring set.
- [x] Attribute sample-narrowed MS evidence only to the best predicted allele
      within the donor allele set.
- [x] Keep genuinely class-only rows with no usable exact or donor-set
      assignment on the stricter ``unrestricted_ms`` path.
- [x] Add regression coverage for class-only donor-set evidence and
      ``peptide_attribution`` provenance.
- [x] Update docs to clarify that sample-bag rows are prediction-deconvolved.
- [x] Bump the patch version.
- [x] Run ``./format.sh``, ``./lint.sh``, and ``./test.sh``.
- [ ] Merge and deploy the PR.

## Review

- ``sample_allele_match`` and ``peptide_attribution`` now share the
  sample-narrowed evidence path.
- Class-only rows with donor allele bags are deconvolved by prediction: only
  the best predicted allele in the bag gets ``sample_allele_ms`` support.
- Rows with no usable exact or donor-set allele assignment remain on the
  stricter ``unrestricted_ms`` path.
- Documentation now states that donor bags are candidate sets, not direct HLA
  assignments.
- Verification passed: ``./format.sh``, ``./lint.sh``, and ``./test.sh``
  (306 tests).

---

# PR - Use Live Hitlist MS Evidence For Panel Construction (2026-05-06)

## Goal

Stop treating packaged gene-level MS counts as the source of truth for automatic
CTA panel construction. Panel selection should rank by non-MS signals first, then
compute MS safety and CTA-exclusive support for candidate CTAs from the current
hitlist observations index. The packaged CSV should not need bundled
``ms_*``/``ms_cta_exclusive_*`` count columns for panel behavior.

## Plan

- [x] Make automatic panel selection ignore bundled MS count/confidence columns
      and use live hitlist-derived MS support for candidate filtering/ranking.
- [x] Add a clear error when automatic panel construction needs public MS
      evidence but the hitlist index cannot be built or read.
- [x] Remove packaged MS count columns from ``cancer-testis-antigens.csv`` and
      adjust evidence helpers/tests to treat those counts as runtime-derived.
- [x] Preserve explicit CTA requests and existing peptide-level pMHC evidence
      behavior.
- [x] Update docs to say hitlist observations are a panel dependency, not an
      optional freshness improvement.
- [x] Bump the patch version.
- [x] Run ``./format.sh``, ``./lint.sh``, and ``./test.sh``.
- [ ] Open, merge, and deploy the PR.

## Review

- Automatic panel construction now ranks initial candidates by non-MS columns
  such as HPA tumor prevalence, then recomputes current hitlist-derived public
  MS support, live MS restriction confidence, and vital healthy-MS vetoes for
  each candidate batch before pMHC scoring.
- ``tsarina/data/cancer-testis-antigens.csv`` no longer bundles runtime
  ``ms_*count*`` columns; tests assert packaged evidence remains count-free.
- Removed stale MS-count rank behavior: asking to rank by removed packaged
  ``ms_*count*`` columns raises a clear error instead of silently falling back
  to alphabetical ordering.
- Explicit CTA requests still bypass automatic CTA-family gates, while the
  per-peptide pMHC evidence path continues to use current hitlist observations.
- Verification passed: ``./format.sh``, ``./lint.sh``, and ``./test.sh``
  (304 tests).

---

# PR - Fix XAGE CTA-Exclusive MS Counts (2026-05-06)

## Goal

Fix stale bundled XAGE1A/XAGE1B MS evidence counts so shared XAGE peptides that
are absent from non-CTA proteins count as CTA-exclusive evidence for both
paralogs and for the grouped XAGE target. Also regenerate stale
``restriction_confidence`` values so PAGE2/PAGE5 and other genes match the
current synthesis algorithm. Track this as
https://github.com/pirl-unc/tsarina/issues/58 and
https://github.com/pirl-unc/tsarina/issues/59.

## Plan

- [x] Update bundled XAGE1A/XAGE1B rows with the recomputed cancer-only MS and
      CTA-exclusive cancer-MS counts.
- [x] Regenerate bundled ``restriction_confidence`` values from the current
      synthesis algorithm.
- [x] Add a regression test locking the expected XAGE shared-peptide count.
- [x] Extend the tier consistency test to check ``restriction_confidence``.
- [x] Add compact HPA cancer RNA/IHC prevalence tables for the full bundled CTA
      candidate universe.
- [x] Add a regeneration script that streams HPA RNA sample data and subsets it
      to whatever CTA CSV is present, so larger future CTA sets are handled.
- [x] Use bundled tumor RNA prevalence features in automatic panel ranking while
      keeping the previous MS-count rank available via ``--cta-rank-by``.
- [x] Bump the patch version.
- [x] Run ``./format.sh``, ``./lint.sh``, and ``./test.sh``.
- [ ] Open, merge, and deploy the PR.

## Review

- Recomputed XAGE1A/XAGE1B packaged MS evidence so their four shared
  cancer-MS peptides count as CTA-exclusive for both paralogs.
- Regenerated stale ``restriction_confidence`` values and extended the
  consistency test to catch future drift.
- Added bundled HPA cancer RNA/IHC prevalence summaries for all 358 current CTA
  candidates plus ``scripts/regenerate_hpa_cancer_prevalence.py`` for future
  larger candidate sets.
- Automatic panel ranking now defaults to HPA tumor RNA prevalence breadth /
  sample prevalence, with MS-supported CTAs sorted before zero-MS candidates and
  the previous MS-first ranking still available via
  ``--cta-rank-by ms_cta_exclusive_cancer_peptide_count``.
- Real default-panel smoke run selected 25 targets including
  ``XAGE1A/XAGE1B``, ``PAGE2/PAGE2B``, and ``PAGE5``.
- Verification passed: ``./format.sh``, ``./lint.sh``, and ``./test.sh``
  (301 passed, 1 skipped locally because the optional sibling ``mhcflurry``
  checkout is not importable).

---

# PR - Avoid Default MHCflurry Allele Chunking (2026-05-05)

## Goal

Fix the interactive ``tsarina panel`` slow path where tqdm progress bars chunk
MHCflurry scoring by allele and repeat the allele-independent processing-model
pass many times. Track this as https://github.com/pirl-unc/tsarina/issues/56.

## Plan

- [x] Default MHCflurry progress-bar scoring to one full allele batch, while
      preserving explicit ``--score-chunk-size`` overrides.
- [x] Update CLI/docs/tests to explain that MHCflurry chunking is opt-in
      because it can repeat processing work.
- [x] Bump the patch version.
- [x] Run ``./format.sh``, ``./lint.sh``, and ``./test.sh``.
- [ ] Open, merge, and deploy the PR.

## Review

- MHCflurry progress-bar scoring now defaults to one full allele batch, avoiding
  repeated processing-model passes across allele chunks.
- Explicit ``--score-chunk-size`` still works for MHCflurry when a user chooses
  to trade speed for smaller scoring chunks.
- CLI/docs now describe the MHCflurry chunking tradeoff.
- Verification passed: ``./format.sh``, ``./lint.sh``, and ``./test.sh``
  (297 tests).

---

# PR - Display CTAG1A/CTAG1B As The Canonical Group Label (2026-05-05)

## Goal

Use the gene-pair display label ``CTAG1A/CTAG1B`` for the grouped CTAG1 target
instead of ``NY-ESO-1``, while continuing to accept ``NY-ESO-1`` as a clinical
alias on input.

## Plan

- [x] Change the canonical CTAG1 group label in the spanning selector and CLI
      defaults.
- [x] Update docs and tests so output names use ``CTAG1A/CTAG1B`` and aliases
      still normalize correctly.
- [x] Bump the patch version.
- [x] Run ``./format.sh``, ``./lint.sh``, and ``./test.sh``.
- [ ] Open, merge, and deploy the PR.

## Review

- Canonical CTAG1 output now uses ``CTAG1A/CTAG1B`` in the default allowlist,
  automatic panel rows, CLI defaults, and docs.
- ``NY-ESO-1`` remains an accepted input alias and still expands to both
  ``CTAG1A`` and ``CTAG1B`` for peptide enumeration.
- Verification passed: ``./format.sh``, ``./lint.sh``, and ``./test.sh``
  (296 tests).

---

# PR — Group CTA Targets With Identical Peptide Sets (2026-05-05)

## Goal

Group paralogous CTA target labels when their enumerated CTA-exclusive peptide
sets are identical, so families such as XAGE/PAGE can appear as combined panel
target names when they genuinely collapse to the same peptide set. Keep
``NY-ESO-1`` as the combined CTAG1A/CTAG1B display target and make progress
messages distinguish requested target labels from underlying Ensembl genes.

## Plan

- [x] Add peptide-set grouping metadata during CTA peptide enumeration without
      losing the original source target labels.
- [x] Teach automatic backfill and output grouping to count peptide-identical
      CTA target groups before falling back to final pMHC-signature grouping.
- [x] Preserve combined names and source members in long output metadata and
      summaries.
- [x] Clarify CLI/docs progress semantics for target labels versus Ensembl gene
      expansion.
- [x] Add regression tests for XAGE/PAGE-style peptide-set grouping and for the
      expanded-gene progress message.
- [x] Bump the patch version and run ``./format.sh``, ``./lint.sh``, and
      ``./test.sh``.

## Review

- Added peptide-set grouping metadata during CTA peptide enumeration. Source CTA
  labels remain intact until output grouping, so automatic backfill order still
  follows the ranked candidate list.
- Peptide-identical paralog targets now appear under combined names such as
  ``XAGE1A/XAGE1B`` when their enumerated peptide sets match exactly.
- Added ``--no-group-identical-cta-peptide-sets`` for users who want separate
  paralog rows despite identical peptide sets.
- Progress now reports when target labels expand to multiple Ensembl genes, so
  ``NY-ESO-1`` can explain a 25-target / 26-gene enumeration step.
- Verification passed: ``./format.sh``, ``./lint.sh``, and ``./test.sh``
  (296 tests).

---

# PR — Group Redundant CTA Panels And NetMHCpan Affinity Annotation (2026-05-05)

## Goal

Collapse CTA targets that yield identical selected pMHC panels so paralogous
families do not consume multiple automatic panel slots, clarify broad
``unrestricted_ms`` evidence, and optionally annotate every selected pMHC with
NetMHCpan binding affinity and affinity percentile rank.

## Plan

- [x] Group non-empty CTAs by their final selected pMHC signature, preserving
      member-gene provenance in output metadata and long tables.
- [x] Make automatic ``cta_count`` count grouped non-empty CTA targets, so
      duplicate peptide panels backfill with additional distinct targets.
- [x] Add optional NetMHCpan BA scoring for selected pMHCs, emitting
      NetMHCpan affinity nM and affinity percentile columns without changing
      the default predictor path.
- [x] Preserve affinity percentile output from topiary/mhctools predictors
      when available.
- [x] Update CLI, docs, and tests for grouped CTAs, ``unrestricted_ms`` meaning,
      and the NetMHCpan annotation option.
- [x] Bump the patch version and run ``./format.sh``, ``./lint.sh``, and
      ``./test.sh``.

## Review

- Added default grouping for CTAs that have identical final selected pMHC
  signatures; grouped outputs preserve source labels in ``cta_members`` and
  ``cta_groups`` metadata.
- Automatic panel backfill now counts grouped non-empty CTA targets, so duplicate
  paralog panels do not consume multiple requested slots.
- Added ``--netmhcpan-affinity`` / ``annotate_netmhcpan_affinity=True`` to add
  NetMHCpan BA nM and affinity percentile rank columns for selected pMHCs.
- Preserved affinity percentile ranks from topiary/mhctools predictors in the
  shared scoring API; MHCflurry still avoids affinity percentile calibration and
  returns that column empty.
- Clarified that ``unrestricted_ms`` means class-I peptide-level MS evidence
  with no usable allele assignment; selected HLAs are prediction-assigned under
  the stricter unrestricted-MS cutoff.
- Verification passed: ``./format.sh``, ``./lint.sh``, and ``./test.sh``
  (294 tests).

---

# PR — CTA Filter Threshold And Column Rename (2026-05-05)

## Goal

Relax the no-protein / uncertain-protein CTA RNA gate from 99% to 98% so
borderline testis-restricted candidates like `XAGE1B` pass, and rename the
bundled evidence inclusion column from `filtered` to `passes_filters`.

## Plan

- [x] Change the adaptive HPA RNA threshold for `Uncertain` / missing protein
      evidence from `0.99` to `0.98`.
- [x] Rename the bundled CTA evidence column to `passes_filters`, while keeping
      existing `CTA_filtered_*` helper names as compatibility aliases.
- [x] Update code paths that filter CTA evidence to use a shared
      `passes_filters` mask and tolerate older CSVs with `filtered`.
- [x] Regenerate or mechanically update the bundled CTA evidence table so
      `XAGE1B` passes and counts move from 257/278 to 258/279.
- [x] Update docs and tests for the new threshold, column name, and counts.
- [x] Bump the package patch version and run `./format.sh`, `./lint.sh`, and
      `./test.sh`.

## Review

- Lowered the missing/uncertain-protein adaptive RNA threshold to `0.98`,
  adding `XAGE1B` to the expressed CTA set.
- Renamed the bundled evidence column to `passes_filters`; internal filtering
  now uses a shared mask helper that still accepts legacy `filtered` inputs.
- Added `rna_98_pct_filter` to the bundled table and updated docs/tests/counts
  from 257/278 to 258/279.
- Verification passed: `./format.sh`, `./lint.sh`, and `./test.sh`
  (`290 passed`).

---

# PR — CTA-Exclusive MS Evidence Counts (2026-05-05)

## Goal

Regenerate the bundled CTA evidence table with CTA-exclusive MS peptide counts,
rank default panel targets by the CTA-exclusive cancer-MS count, and document
why PAGE/XAGE peptides pass or fail the strict CTA-exclusive panel path.

## Plan

- [x] Add CTA-exclusive MS count columns to the gene-level CTA evidence
      aggregation while preserving the existing all-CTA MS count columns.
- [x] Recompute `tsarina/data/cancer-testis-antigens.csv` so the bundled table
      matches the new PR behavior.
- [x] Use CTA-exclusive cancer-MS evidence as the default automatic panel rank
      column.
- [x] Compute estimated population coverage by summing covered allele
      frequencies within each HLA locus before converting to carrier
      probabilities.
- [x] Add regression tests for the split between all-CTA MS evidence and
      CTA-exclusive MS evidence.
- [x] Audit PAGE2, PAGE5, and XAGE1A selected/raw peptides against the protein
      catalogue to explain cross-protein occurrence and CTA-exclusivity.
- [x] Bump the package patch version and run `./format.sh`, `./lint.sh`, and
      `./test.sh`.

## Review

- Added CTA-exclusive MS count columns to the bundled CTA evidence table while
  preserving broad all-CTA MS count columns.
- Regenerated `tsarina/data/cancer-testis-antigens.csv`: PAGE2 now has 9
  CTA-exclusive cancer-MS peptides, PAGE5 has 6, and XAGE1A has 0 because its
  raw MS peptides also occur in XAGE1B, which is outside the clean CTA
  partition.
- Automatic CTA ranking now defaults to
  `ms_cta_exclusive_cancer_peptide_count`; explicit `--cta-rank-by
  ms_cancer_peptide_count` still requests the older broad count.
- Estimated CTA population coverage now sums covered allele frequencies within
  each HLA locus, converts each locus total to carrier probability, then
  combines HLA-A/B/C loci.
- Documented that the cached public-MS path uses hitlist observations built
  from registered IEDB/CEDAR plus hitlist supplementary MS rows, not bulk
  proteomics or line-expression indexes.
- Verification passed: `./format.sh`, `./lint.sh`, and `./test.sh`
  (288 tests).

# PR — Panel MAGE-Family Safety Gate And Summary Sorting (2026-05-05)

## Goal

Make the default CTA panel selection safer around the MAGE family and make
coverage summaries easier to scan by sorting CTA rows by selected peptide yield.

## Plan

- [x] Add an automatic MAGE-family safety gate that allows `MAGEA4` by default
      but excludes other `MAGE*` CTAs unless they are explicitly requested or
      allowlisted.
- [x] Expose a CLI switch to disable the automatic MAGE-family gate for users
      who intentionally want broader MAGE-family exploration.
- [x] Sort "Expected Population Coverage Per CTA" rows by selected peptide
      count, then HLA hits, then estimated coverage, while preserving zero-hit
      rows so failed candidates are visible.
- [x] Add numeric frequency support for every default-panel allele so selected
      HLA hits do not report artificial zero estimated coverage.
- [x] Add a frequency audit layer that keeps regional proxy frequencies and
      published global averages on the same allele-frequency scale, records
      source/proxy/resolution provenance, and verifies every default-panel
      allele has a global average plus a coverage frequency.
- [x] Pin clinical allowlisted CTAs into automatic panels and hide downstream
      empty CTAs from default automatic panel output unless explicitly requested.
- [x] Ensure sample-genotype MS evidence assigns HLA specificity by
      best-of-haplotype prediction unless the evidence is monoallelic.
- [x] Backfill automatic CTA panels from lower-ranked candidates until the
      requested count of downstream non-empty CTAs is reached when possible.
- [x] Split monoallelic MS pMHC support from sample/deconvolved MS support in
      per-CTA coverage summaries.
- [x] Document why selected CTAs may have zero peptides after downstream
      peptide/exclusivity/MS/prediction gates.
- [x] Bump the package patch version and run `./format.sh`, `./lint.sh`, and
      `./test.sh`.

## Review

- Default automatic CTA selection now excludes `MAGE*` targets other than
  `MAGEA4` unless they are explicitly selected or allowlisted.
- Added `--allow-non-magea4-mage-family` for deliberate broader MAGE-family
  exploration.
- "Expected Population Coverage Per CTA" rows now sort by selected peptide
  count, then HLA-hit count, then estimated coverage.
- All `global53_abc` alleles now have numeric frequency support: regional proxy
  rows when available, otherwise published global CIWD fallbacks; `C*04:03`
  uses a specific East Asian AFND proxy row.
- Added an allele-frequency audit table that exposes regional weighted
  frequency, published global average, coverage frequency, coverage source, and
  exact/proxy/qualitative regional support counts. All 53 default-panel alleles
  have published global averages; coverage uses regional weighted frequencies
  for 35 and published global averages for 18 with no numeric regional proxy.
- Automatic panel selection now pins the default clinical allowlist
  (`MAGEA4`, `PRAME`, `NY-ESO-1`) ahead of lower-ranked candidates, and default
  automatic output hides CTAs with no selected pMHCs. Use `--show-empty-ctas`
  to restore the previous audit view; explicit `--ctas` requests are preserved
  even when empty.
- Sample-genotype MS rows now use best-of-haplotype specificity before broad
  exact-restriction assignment, so multi-allelic sample rows that list several
  HLA restrictions only support the best predicted panel allele unless they are
  monoallelic.
- Automatic panel selection now scans lower-ranked candidates in batches to
  backfill downstream-empty CTAs, so `cta_count=25` means up to 25 non-empty
  selected CTA targets when enough candidates pass downstream peptide/MS/HLA
  gates. `--show-empty-ctas` still restores the top-candidate audit view.
- Per-CTA coverage summary rows now split selected pMHC counts into
  monoallelic MS, sample/deconvolved MS, and unrestricted MS columns.
- Verification passed: `./format.sh`, `./lint.sh`, and `./test.sh`
  (286 tests).

# PR — Faster MHCflurry Scoring Without Affinity Percentile Calibration (2026-05-04)

## Goal

Fix MHCflurry failures for alleles such as `HLA-C*15:05` where affinity
prediction works but affinity percentile-rank calibration is missing, and reduce
panel/personalization scoring time by avoiding unnecessary predictor work.

## Plan

- [x] Add a direct MHCflurry scoring path that returns tsarina's existing
      `presentation_score`, `presentation_percentile`, and `affinity_nm`
      columns without requesting unused affinity percentile ranks.
- [x] Batch MHCflurry presentation predictions across peptide/allele pairs
      instead of calling the presentation predictor once per allele.
- [x] Keep non-MHCflurry predictors on the existing topiary/mhctools path.
- [x] Add regression tests for missing affinity percentile calibration and
      batched direct-MHCflurry calls.
- [x] Bump the package patch version.
- [x] Run `./format.sh`, `./lint.sh`, and `./test.sh`.

## Review

- Filed upstream mhctools wrapper issue:
  https://github.com/openvax/mhctools/issues/203.
- Local `HLA-C*15:05` MHCflurry smoke test returns `presentation_score`,
  `presentation_percentile`, and `affinity_nm` without requesting affinity
  percentile ranks.
- Warm local timing on 24 peptide-allele pairs: direct MHCflurry path `0.063s`
  vs old topiary/mhctools MHCflurry path `0.500s`.
- Verification passed: `./format.sh`, `./lint.sh`, and `./test.sh` (260 tests).

## Follow-Up Plan — Audit Global-51 Allele Panel

- [x] Cross-check the existing `global51_abc_ssa` panel against MHCflurry's
      runtime affinity percentile-rank calibration resolver.
- [x] Compare the panel against IEDB/TepiTool global class-I reference-set
      alleles, the IEDB/Paul 38 common A/B threshold panel, and Sarkizova
      HLA-C frequent allotypes.
- [x] Replace uncalibrated or weakly justified add-ons only when a stronger
      calibrated, reference-backed allele is available.
- [x] Record citation/provenance notes in `PANEL_SOURCE_CATEGORIES`.
- [x] Add tests that keep `global51_abc` MHCflurry-compatible and prevent
      future weak local-only complements from entering the default panel.
- [x] Rerun `./format.sh`, `./lint.sh`, and `./test.sh`, then update PR #47.

## Review — Audit Global-51 Allele Panel

- Existing `global51_abc_ssa` works with MHCflurry's runtime percentile-rank
  calibration resolver, including `HLA-A*24:02`.
- Added `global51_abc` as the new default panel. It uses all 27
  IEDB/TepiTool A/B alleles, all 21 Sarkizova frequent HLA-C allotypes, and
  three highest-frequency calibrated
  IEDB/Paul common A/B complements.
- Excluded `HLA-C*15:05` because MHCflurry supports raw affinity and presentation
  prediction for it but does not have affinity percentile-rank calibration.
- Verification passed: `./format.sh`, `./lint.sh`, and `./test.sh` (266 tests).

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
- [x] Open PR linked to tsarina#43; merge when CI is green and deploy.

Review:

- Local validation passed: `./format.sh`, `./lint.sh`, and `./test.sh`
  (`256 passed`). Live selector check confirms `NY-ESO-1` remains selected
  without the allowlist at the new 2.0 nTPM threshold, while MAGEA4 remains
  gated by gene-level healthy-MS heart evidence unless allowlisted.
- Shipped as PR #45 / `tsarina` v1.2.9.

---

## Fix Panel Healthy-MS Gene Veto (tsarina#44)

Goal: stop automatic CTA selection from blacklisting a CTA as a whole based on
healthy-MS evidence from peptides that map to multiple CTA-family genes.

Plan:

- [x] Change the panel safety gate so gene-level healthy-MS aggregates are only
      a suspect list; the actual MS veto requires unique vital healthy-MS rows
      for the CTA symbol.
- [x] Keep the RNA vital-tissue threshold behavior from #43 unchanged.
- [x] Add tests showing shared MAGE-family heart MS does not veto MAGEA4, while
      unique MAGEA1 vital MS still gates MAGEA1.
- [x] Update CLI/docs wording from generic healthy-MS veto to unique healthy-MS
      veto.
- [x] Bump version for the PR.
- [x] Run `./format.sh`, `./lint.sh`, and `./test.sh`.
- [ ] Open PR linked to tsarina#44; merge when CI is green and deploy.

Review:

- Local validation passed: `pytest tests/test_spanning.py tests/test_cli_spanning.py -q`
  (`49 passed`), `./format.sh`, `./lint.sh`, and `./test.sh` (`257 passed`).
  Live selector check confirms `MAGEA4` is selected and `MAGEA1` is excluded
  without the allowlist, matching shared-vs-unique healthy-MS evidence.

---

## Rebalance Global Panel With CTA-MS Supported Alleles (PR #47 follow-up)

Goal: add the strongest missing CTA-MS supported alleles (`HLA-B*15:02`,
`HLA-B*27:05`, `HLA-A*29:02`) while preserving the global coverage rationale
of the default panel.

Plan:

- [x] Audit zero-MS `global51_abc` alleles for weighted regional frequency,
      literature/source rationale, and MHCflurry pseudosequence redundancy.
- [x] Decide whether to replace weak alleles or grow the default panel.
- [x] Update panel constants, docs, and tests to reflect the selected panel.
- [x] Run `./format.sh`, `./lint.sh`, and `./test.sh`.
- [x] Push an update to PR #47.

Review:

- Zero-MS alleles with weak local weighted-frequency evidence include
  `HLA-A*30:02`, `HLA-C*03:02`, `HLA-C*04:03`, `HLA-C*07:04`, and
  `HLA-C*14:03`; however these are either IEDB/TepiTool backbone alleles or
  members of the Sarkizova frequent-HLA-C set. `HLA-C*14:03` is redundant with
  `HLA-C*14:02` by MHCflurry pseudosequence, but removing it would break the
  published frequent-C-set rationale.
- Chose to keep `global51_abc` as the 51-allele reference panel and add
  `global53_abc` as the default, adding `HLA-A*29:02`, `HLA-B*15:02`, and
  `HLA-B*27:05` while dropping `HLA-C*14:03` from the default. `HLA-C*14:03`
  is redundant with `HLA-C*14:02` in MHCflurry's pseudosequence and runtime
  percentile-rank calibration, and `HLA-C*14:02` had the CTA-MS support.
- Local validation passed: targeted panel tests (`19 passed`), `./format.sh`,
  `./lint.sh`, and `./test.sh` (`269 passed`).

---

# PR - Add XAGE2 + Lower CTA No-Protein Threshold to 0.97 + Parameterize Floor (2026-06-04)

## Goal

Surface the XAGE-family CTAs the user asked for, properly. Investigation showed
they were in three states: XAGE3 already expressed; XAGE5 passes but flagged
`never_expressed` (1.1 nTPM); XAGE2 absent from the source universe.

## Plan / Decisions (with the user)

- [x] **Add XAGE2** (`ENSG00000155622`, CTpedia/CT12.2; tsarina#79). Built its
      row from HPA `rna_tissue_consensus.tsv` via the in-repo generators
      (`tiers.enrich_rna_per_tissue` + `tiers.assign_all_axes`); validated the
      deflated-fraction formula reproduces XAGE3 (0.9922) and XAGE5 (1.0) and
      that the `passes_filters`/`never_expressed` rules reproduce all 358 shipped
      rows exactly.
- [x] **Investigated XAGE2's lung 5.3 nTPM**: not contamination — replicates
      across GTEx (2.0), HPA (5.3), FANTOM (12.6). Genuine low-level lung
      expression; XAGE2 carries `safety_flags=lung`.
- [x] **Lowered the no-protein/Uncertain adaptive threshold 0.98 -> 0.97**
      (`HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS`). Flips exactly CT83 (KK-LC-1) and
      PRM3 to pass; XAGE2 (0.977) passes. DPPA3 (0.9688) stays out.
- [x] **Parameterized the `never_expressed` floor** as
      `HPA_EXPRESSION_FLOOR_NTPM = 2.0` and **rescued XAGE5** via
      `MANUALLY_EXPRESSED_CTA` rather than a blanket floor drop (tsarina#78).
- [x] Updated docs/curation.md, tests (counts 258->262, 279->282; new XAGE2 /
      XAGE5 / CT83+PRM3 / parameterization tests), minor version bump 1.3.6 ->
      1.4.0.

## Review

- `CTA_gene_names()` 258 -> 262 (+XAGE2, +CT83, +PRM3, +XAGE5). All four carry
  honest evidence: CT83/PRM3/XAGE2 are SOMATIC-restricted with visible
  somatic signals; XAGE2 keeps `safety_flags=lung`.
- XAGE5 stays `never_expressed=True` in the table (HPA truth) but is returned
  by `CTA_gene_names()` via the explicit rescue set — closes tsarina#78.
- `scripts/add_xage2.py` is the reproducible, validated generator (kept in-repo).
- Verified end-to-end: XAGE2/CT83/PRM3 surface in pirlygenes automatically
  (pure data); XAGE5 needs pirlygenes to delegate to tsarina (filed as a
  pirlygenes issue, not implemented here).
- Verification passed: `./format.sh`, `./lint.sh`, `./test.sh` (321 passed).

---

# PR - Add MAGEB6 (CTA), evaluate #79 dual-corroborated candidates (2026-06-04)

## Goal

Pick up tsarina#79: add the CTA genes "corroborated by both CTexploreR and
CTdatabase" that are missing from the panel.

## Finding

Evaluated all 9 dual-corroborated candidates against the real HPA filter
(RNA reproductive-restriction + the protein-in-somatic rule). They are strong
by *database membership* but only **MAGEB6** passes tsarina's HPA standard:

- MAGEB6 (ENSG00000176746): testis-only 3.0 nTPM, no protein -> PASSES (TESTIS).
- RNF17: perfect RNA (1.0) but HPA protein in blood vessel/heart/kidney -> fails.
- TAF7L: protein in pancreas; ROPN1 0.85 + salivary; NXF2 0.37; CT45A5/A6 brain;
  NLRP4 protein in blood vessel -> all fail somatic.
- DSCR8: HPA 404 (retired / not in HPA).

This is exactly the value tsarina's HPA filter adds over raw DB membership.

## Plan

- [x] Add MAGEB6 via `scripts/add_cta_gene.py` (generalizes add_xage2.py;
      asserts each gene passes before writing; idempotent).
- [x] Tests (counts 262->263 / 282->283; MAGEB6 membership) + version 1.4.2->1.4.3.
- [x] Comment the full evaluation on tsarina#79.

## Review

- `CTA_gene_names()` 262 -> 263 (+MAGEB6, clean TESTIS).
- Documented why the other 8 candidates fail on #79 rather than bulk-adding.
- Verification: `./format.sh`, `./lint.sh`, `./test.sh` (323 passed).

---

# PR - Require hitlist 1.55.2 and migrate stale peptide mappings (2026-09-02)

## Goal

Make tsarina's supported class-II gene queries use hitlist's corrected
length-independent peptide mappings, including for users whose existing
`peptide_mappings.parquet` predates hitlist 1.55.2.

## Plan

- [x] Raise the runtime dependency floor from `hitlist>=1.45.0` to
      `hitlist>=1.55.2`.
- [x] Add a cheap, persistent compatibility check for existing mapping
      sidecars. Since hitlist 1.55.2 does not stamp its package/builder version
      into `peptide_mappings_meta.json`, verify the behavior using small,
      deterministic class-II and length-7 observation samples, then cache the
      result against the sidecar's size and nanosecond mtime.
- [x] If observations exist but mappings are missing or fail the compatibility
      probe, rebuild only `peptide_mappings.parquet` with
      `build_peptide_mappings(force=True)` instead of rescanning all evidence.
      Preserve `build_observations(force=True)` for an explicit full rebuild.
- [x] Add regression tests for current, stale, missing, force-built, and
      already-verified caches, including a guard that mapping multi-rows do not
      multiply observation rows.
- [x] Document the automatic one-time migration and its progress message.
- [x] Bump tsarina's patch version for the PR.
- [x] File the missing mapping-builder-version/cache-invalidation contract on
      hitlist and link it from the tsarina PR.
- [x] Run targeted tests, `./format.sh`, `./lint.sh`, and `./test.sh`; review the
      diff for minimality and record results below.
- [ ] Open the PR, merge it after checks pass, then deploy the merged release to
      PyPI from a clean `main` using `./deploy.sh`.

## Review

- Raised the dependency floor to hitlist 1.55.2 and the tsarina patch version
  to 1.24.2.
- Existing current artifacts are behavior-probed once using bounded class-II
  and length-7 samples, then recorded against the mapping parquet's size and
  nanosecond mtime. Replacing the parquet invalidates the marker. The real
  5,862,627-row artifact verified in 2.066 seconds; the fingerprinted repeat
  check took less than 1 millisecond.
- Missing or legacy mappings now call `build_peptide_mappings(force=True)`;
  fresh or explicitly forced observation builds retain the full builder path.
  Gene-filtered `load_ms_evidence` calls always pass through this validation.
- Added regression coverage for all migration branches and confirmed that
  multi-mapping annotation preserves observation-row counts. Targeted suite:
  32 passed. Required gates: `./format.sh`, `./lint.sh`, and `./test.sh`
  (429 passed, 6 existing pandas warnings).
- Filed the upstream root-cause contract as
  https://github.com/pirl-unc/hitlist/issues/404.

## Task: Retire hitlist workarounds superseded by hitlist 1.59.1

hitlist moved 1.55.2 -> 1.59.1 (23 releases). Four pieces of tsarina now
duplicate or fight upstream behavior that hitlist owns, and one dependency
floor no longer describes what tsarina needs.

Verified upstream facts driving this:

- hitlist#44 (`allele_to_serotype` preferred Bw4/Bw6 over the locus serotype)
  closed 2026-04-18. `allele_to_serotype("HLA-A*23:01")` now returns
  `HLA-A23`, `allele_to_all_serotypes` returns `("HLA-A23", "HLA-Bw4")`, the
  observations parquet stores both `serotype` and `serotypes`, and
  `load_observations(serotype=...)` filters on set membership in `serotypes`.
- `allele_resolution` is a stored column. hitlist 1.55.7 made the scanner
  recompute the whole MHC annotation when a class-only row is promoted to a
  donor set, so a promoted row now stores `donor_set` where tsarina's
  per-row recompute of the joined restriction string disagrees with the
  pre-1.55.7 stored value.
- hitlist#404 (the missing mapping-artifact contract that tsarina#147 worked
  around with a behavior probe) is fixed: `_MAPPING_ARTIFACT_VERSION = 2`
  plus a full `contract` block in `peptide_mappings_meta.json`. hitlist#429
  additionally fingerprints the curation YAMLs into
  `observations_meta.json`, and `build_observations(force=False)` /
  `build_peptide_mappings(force=False)` short-circuit when both artifacts
  are valid.
- `restriction_evidence` (hitlist#415) is a study-level claim about how a
  named restriction was established (`experimental` / `monoallelic` /
  `predicted` / `unknown`). It is NOT the same axis as tsarina's
  panel-relative evidence tiers, so `_build_evidence_stats` stays.
  `MHC_ALLELE_PROVENANCE_VALUES` is the new authoritative provenance
  vocabulary that tsarina's hardcoded subset should be checked against.

### Plan

- [x] Replace `_filter_by_serotype`'s mhcgnomes expansion with hitlist's
      `serotypes` column membership, keeping the accepted query spellings
      (`A2`, `HLA-A2`, `Bw4`) and the current match results.
- [x] Make `_apply_min_resolution` read the stored `allele_resolution`
      column instead of reclassifying `mhc_restriction` per row, and let
      `--min-resolution` reach `donor_set`.
- [x] Delete the `.tsarina-peptide-mappings.json` marker and the
      behavior probe; delegate both artifacts' staleness to
      `build_observations(force=...)`, keeping hitlist's stdout off
      tsarina's stdout and tolerating an unregistered IEDB/CEDAR source.
- [x] Drop the stale provenance narrative in `spanning.py` and add drift
      guards that fail if hitlist's provenance / resolution vocabularies
      stop containing the values tsarina keys on.
- [x] Raise the `hitlist` floor to the version that actually provides the
      above and refresh the comment.
- [x] Extend `tests/fixtures/hitlist_mini/*.parquet` with the columns the
      current builder emits (`restriction_evidence`, `gene_biotype`) so the
      one non-mocked integration test can catch schema drift.
- [x] Bump the tsarina version.
- [x] File the upstream gaps this PR has to work around (public serotype
      query normalizer, public cache-validity accessor) on hitlist and link
      them from the PR.
- [x] Run `./format.sh`, `./lint.sh`, `./test.sh`; record results below.
- [ ] Open the PR, merge after checks pass, then `./deploy.sh` from clean
      `main`.

### Review

- `--serotype` now reads hitlist's `serotypes` column instead of expanding the
  query through mhcgnomes, deleting the hitlist#44 workaround. NOTE: the two
  behavior changes claimed here were wrong and are corrected in the follow-up
  task below — public epitopes already worked before this PR, and matching
  donor sets was a mistake, reverted.
- `--min-resolution` reads the stored `allele_resolution` instead of
  reclassifying the restriction string, and `donor_set` became reachable (it
  ranks between `four_digit` and `two_digit` and was previously unofferable).
  Both filters now fail with an actionable message on an index too old to carry
  the column, rather than answering from a re-derivation.
- `ensure_index_built` hands freshness to `build_observations(force=...)`,
  which fingerprints the curation YAMLs and both artifact contracts. That
  retires the `.tsarina-peptide-mappings.json` marker and the two-probe
  behavior test added in #147 as a workaround for hitlist#404, now fixed
  upstream. Existence-only gating was letting a pre-#412 artifact answer every
  query: on this machine that artifact still carried 942 purified-MHC /
  half-life rows as MS evidence that hitlist#423 reclassifies as binding
  assays, plus retired sample attributions for four corrected studies.
- hitlist prints its build report to stdout, which is `tsarina hits`' data
  stream. A known build now streams to stderr; a validation run captures the
  report and replays it only if an artifact actually changed, so a current
  cache stays silent. Filed pirl-unc/hitlist#448 for the public, quiet
  validity predicate that would let tsarina announce a stale rebuild before
  spending ten minutes on it, and pirl-unc/hitlist#449 for the serotype query
  normalizer that `_serotype_query` currently copies.
- An index copied in without its IEDB/CEDAR exports cannot be fingerprinted, so
  it is used as found. `tests/test_hitlist_integration.py` exercises exactly
  that path against the committed fixture.
- `_SAMPLE_NARROWED_PROVENANCES` keeps its two values — `restriction_evidence`
  (hitlist#415) is a study-level claim about how a restriction was established,
  not the panel-relative axis `_build_evidence_stats` computes, so the tier
  derivation stays. `tests/test_hitlist_vocabularies.py` now fails if hitlist's
  provenance or resolution vocabulary moves under either literal.
- `scripts/regenerate_hitlist_mini_fixture.py` makes the fixture reproducible
  for the first time (its slice was previously recoverable only by
  reverse-engineering the peptide list) and `--check` reports drift. The
  refreshed slice gains `restriction_evidence`, `gene_biotype`, `cell_type` and
  the scanner's MHC identity block, and drops two columns hitlist derives at
  load time rather than storing. Two new tests assert the fixture carries what
  a current build writes and that its stored annotations still match the
  installed hitlist.
- Gates: `./format.sh`, `./lint.sh` clean; `./test.sh` 442 passed (was 429),
  6 pre-existing pandas warnings.

## Task: One coding-gene universe for all of tsarina

`viral.py` disagrees with itself about what "human" means, and the disagreement
is not local to `viral.py` — it is four independent restatements of "which
Ensembl biotypes count as coding", at two different levels.

hitlist 1.55.8 widened `proteome_kmer_set`'s default from `protein_coding` to
`ENSEMBL_CODING_GENE_BIOTYPES` (adds IG_V/D/J/C and TR_V/D/J/C germline
segments). tsarina picked that up silently in one code path and not the others.

### The four sites

Gene level:
- `partition.py:76` — `g.biotype == "protein_coding"`, builds the CTA / non-CTA
  partition. 20,089 genes against hitlist's 20,500.

Transcript level:
- `peptides.py:185` — canonical transcript for CTA peptide enumeration
- `peptides.py:275` — non-CTA peptide enumeration, the specificity screen
- `qc.py:75` — longest protein length for fragment-gene-model QC

The two levels must move together. Ensembl gives an IG_V gene's transcripts the
biotype `IG_V_gene`, not `protein_coding`, so widening the gene universe alone
would add 411 genes that contribute zero peptides — an inert change. Verified
those transcripts carry real translations (IGKV4-1 121 aa, TRGV11 103 aa), so
they are genuinely presentable.

### Why widen rather than narrow

Germline IG/TR segments are expressed self peptides, which is why hitlist
changed its default. For a therapeutic target screen the consequences both run
the safe way: a larger self set drops more viral peptides as human, and treating
IG/TR as non-CTA makes a CTA peptide sharing a sequence with an IG/TR segment
fail specificity rather than pass it.

### Plan

- [x] One definition, sourced from hitlist so it cannot drift again:
      `CODING_GENE_BIOTYPES` plus `is_coding_gene` / `is_coding_transcript`
      predicates in `gene_sets.py`, tsarina's existing gene-universe authority
      layer.
- [x] Route all four sites through the predicates.
- [x] Measure: genes gained by the partition, whether any IG/TR gene lands in
      the CTA set, and how many CTA peptides newly fail specificity.
- [x] Drift guard: tsarina's set must equal hitlist's, and must include the
      IG/TR biotypes by name.
- [x] Version bump, three gates, PR.

### Review

- `gene_sets.CODING_GENE_BIOTYPES` aliases hitlist's
  `ENSEMBL_CODING_GENE_BIOTYPES`, with `is_coding_gene` / `is_coding_transcript`
  predicates. All four sites route through them, and no
  `biotype == "protein_coding"` comparison remains anywhere in the package.
- The transcript-level half was the load-bearing part. Ensembl gives an IG_V
  gene's transcripts the biotype `IG_V_gene`, and those transcripts carry real
  translations (IGKV4-1 121 aa, TRGV11 103 aa), so widening only the gene
  universe would have added 411 genes contributing zero peptides.
- Measured on Ensembl 112: 411 IG/TR coding genes, all 411 landing in `non_cta`
  and **none** in `cta`, taking the partition to 293 CTA / 20,198 non-CTA. Those
  segments contribute 92,780 distinct 8-11mers.
- Blast radius on CTA specificity is **one peptide**: `LEGPLRLS` in CTAGE1
  (ENST00000391403, position 510) also occurs in an IG/TR segment, so it now
  fails CTA-exclusivity instead of passing it. It is an 8-mer -- the shortest
  length enumerated -- which is what a chance collision looks like rather than
  real homology. 417,868 of 417,869 CTA peptides are unaffected.
- Nothing changes for `human_exclusive_viral_peptides`, which already used the
  wide set; what changes is that `cancer_specific_viral_peptides` now partitions
  the same universe, so a viral peptide matching an IG/TR segment is attributable
  to non-CTA instead of falling into neither bucket.
- CTA protein-length QC is untouched: only IG/TR genes have IG/TR transcripts,
  and no CTA is one.
- Gates: `./format.sh`, `./lint.sh` clean; `./test.sh` 445 passed (was 443).

## Task: Correct the serotype filter (follow-up to #148)

#148 replaced the mhcgnomes serotype expansion with membership in hitlist's
`serotypes` column. Two things were wrong with it, found by diffing the old
and new filters over all 980 distinct human restrictions in the index rather
than reasoning about them.

1. **The delta was not what the PR claimed.** Public epitopes were already
   queryable before #148 — the old code expanded `Bw4` into its member alleles,
   so `--serotype Bw4` matched A*23:01 and A*24:02 all along. What actually
   changed was that donor sets started matching, taking `--serotype A2` from 14
   to 272 distinct restrictions. This is candidate membership, consistent
   with `--allele`; it does not credit a donor bag to one presenter before
   deconvolution. Resolution thresholds control which candidate sets to keep.
2. **Lowercase regressed.** mhcgnomes parses case-insensitively, so the old
   filter matched molecular rows for `--serotype bw4` and `hla-a24`. The
   hand-rolled `HLA-` prefix rule that replaced it did not, and its docstring
   claimed otherwise.

### Plan

- [x] Put one cached mhcgnomes parse in `tsarina/mhc.py` that a caller can
      pin to an expected reading (`serotype` vs `allele`), and derive the
      serotype comparison key from it so query and stored token agree by
      construction.
- [x] ~~Exclude donor sets from `--serotype`~~ — reverted before merge. Donor
      sets match on membership, exactly as `--allele` has always matched a
      semicolon-joined restriction. Narrowing to single-molecule restrictions is
      `--min-resolution`'s existing job; a second, differently-drawn line inside
      `--serotype` would have made the two filters disagree about the same row.
- [x] Fail loudly on a serotype query that cannot be read, instead of
      returning every row.
- [x] Re-diff old vs new over the full restriction vocabulary and require the
      only differences to be improvements.
- [x] Correct the overstated claims in the #148 review notes below.
- [x] Bump the version, run the three gates.

### Review

- `tsarina/mhc.py` now owns one `parse_mhc(value, expect=...)`, LRU-cached on
  `(value, expect)`, using mhcgnomes' `required_result_types` so a stated
  expectation — a CLI flag, or a curated paper record that reports a
  serological typing rather than a molecule — decides how an ambiguous token is
  read. `serotype_key` reduces both sides of the comparison through it, so case,
  the `HLA-` prefix, and split serotypes stop being tsarina's problem.
- Three legacy curated serotype names (`DR1B`, `DR3A`, `DR7A`) remain
  queryable through an explicit exception set when mhcgnomes cannot parse
  them. Every other token must parse as a serotype, so compact molecular
  aliases and unknown serotype-shaped labels are rejected.
- Diffed old vs new across all 980 distinct human restrictions for 16 queries.
  Nothing is lost on any query. The only additions: `A*24:03` for `A24` and
  `A*02:03` for `A2` (both are split-serotype members whose broad parent
  mhcgnomes' own member list omits — hitlist's reverse map adds it, and the
  allele genuinely belongs to the parent), `A*24:03` for the split query
  `A2403` itself, and serological rows for lowercase queries, which the old
  filter matched for molecular rows only.
- On donor sets, the first version of this branch excluded them and that was
  wrong. `--allele` matches a donor set containing the queried allele -- by
  design, with a committed test -- so excluding them from `--serotype` made two
  filters answer differently about one row, which is the split-brain the
  coding-universe task had just removed. It also treated the largest rung of the
  corpus as an edge case: on human class I, 1,411,961 rows name one allele and
  1,431,499 give only a candidate set, donor sets being 902,990 of them.
  `--serotype` is therefore a plain membership test, and the existing
  `--min-resolution` is how a caller asks for restrictions that name one
  molecule. Verified composing on the fixture: no filter 16 pMHC rows,
  `--serotype A2` 8, `--serotype A2 --min-resolution four_digit` 3 with no donor
  sets -- the same 3 the exclusion produced, reached with the flag that already
  existed.
- Gates: `./format.sh`, `./lint.sh` clean; `./test.sh` 452 passed.


## Task: Address PR #149 review findings (#153, #154)

### Specification

- Preserve the existing local donor-set membership correction and include it
  in PR #149. `--serotype` and `--allele` must agree on candidate membership;
  `--min-resolution donor_set` keeps donor rows, while `four_digit` removes
  them before the serotype filter runs.
- Replace the serotype-shape regex fallback with an explicit exception set
  for the three previously supported curated names DR1B, DR3A, and DR7A.
  All other tokens must parse as mhcgnomes Serotype results. Preserve case,
  whitespace, and optional HLA-prefix handling for those exceptions.
- Verify both molecular spellings and unknown serotype-shaped queries raise
  even in mixed valid/invalid requests, while valid split and public-epitope
  queries retain their result sets. Do not change allele-filter semantics.
- Keep the PR's existing 1.25.4 patch bump if it remains unreleased. Rewrite
  the PR description around final behavior and link both closing issues.
- Run format, lint, full tests and corpus comparison; merge only the checked
  commit, deploy with `./deploy.sh` from clean main, and verify PyPI artifacts.
- Inspect related open issues after release to identify the next work block.

### Plan

- [x] Re-read current code, outstanding changes, issues, and deployment script.
- [x] Create feature branch and record implementation/verification plan.
- [x] Add regression tests and demonstrate the validation failure (6 failures).
- [x] Narrow the fallback, retain donor-set correction, and update docs/lessons.
- [x] Run full checks and compare the real restriction vocabulary.
- [ ] Commit, update/push PR #149, and wait for CI.
- [ ] Merge and deploy 1.25.4 from clean main; verify publication.
- [ ] Record results and identify the next dependency/urgency work block.

### Review before shipping

- Both review findings are fixed: donor sets match by serotype membership,
  and invalid compact molecular aliases no longer bypass serotype validation.
- The fallback now accepts only the three legacy curated exceptions; tests
  preserve their case/prefix behavior and reject unknown lookalike names.
- The new validation coverage failed in six cases before the fix. All 39
  focused tests pass after it, including donor_set/four_digit composition.
- Corpus comparison: all 142 stored human serotype queries preserve main's
  result sets across 980 distinct restrictions and match lowercase variants.
  All 171 installed HLA serotype table names preserve case/prefix keys.
- Required checks passed: `./format.sh` (unchanged), `./lint.sh`, and
  `./test.sh` (463 passed, 6 warnings).
- Version 1.25.4 is already bumped in this PR; PyPI currently has 1.25.3.
- Final merge/deployment results will be recorded on PR #149 after shipping.

---

## Post-merge /code-review follow-ups on #149 (1.25.5)

`/code-review` on the merged #149 diff (mhc.py, cli_hits.py) surfaced 11
findings, filed nowhere (no issue tracker session open for this) but worked
directly since none required a design call. All 11 fixed on
`fix/code-review-serotype-followups`.

### Plan

- [x] Read the review findings against current `tsarina/mhc.py` /
      `tsarina/cli_hits.py`, confirm each still reproduces.
- [x] Branch off main (`fix/code-review-serotype-followups`).
- [x] Fix all 11: 4 real bugs (uncaught ValueError, KeyError/AttributeError
      footguns in `parse_mhc`, silent cross-species reparsing), 2 efficiency
      (late `--serotype` validation, row-wise serotype match), 2
      simplification (stringly-typed class lookup, dead `_parse_hla`
      passthrough), 1 reuse (duplicated HLA-prefix stripping), 1 docs
      (`serotype_key` docstring overclaim).
- [x] Add/extend regression tests for each fix.
- [x] Run `./format.sh`, `./lint.sh`, `./test.sh`.
- [ ] Commit, push, open PR, merge, deploy, verify PyPI.

### Review

The two most severe findings were confirmed by the reviewer via direct
execution, and stayed reproduced going into the fix:

- **Uncaught crash.** `_filter_by_serotype`'s `ValueError` on an unreadable
  `--serotype` token was never caught between `handle()` and `main()`.
  Fixed two ways: `--serotype` now validates at argparse parse time via a
  new `_parse_serotypes` type= callable (mirrors `--lengths`'
  `_parse_lengths`), so a typo fails before gene resolution or any index
  load runs; and `handle()` also wraps the filter calls in `try/except
  ValueError` for any caller that builds an `argparse.Namespace` directly
  and bypasses the parser.
- **Silent cross-species reparsing.** `normalize_mhc_restriction` lost its
  old hand-rolled `HLA-` force-prefix in #149 (correctly — that prefix
  trick separately regressed lowercase queries) but nothing replaced its
  species-scoping effect, so a non-human token like `SLA1*01:01` was
  silently reparsed and reformatted as a different species' canonical form
  (`SLA-1*01:01`) instead of passed through unchanged. Fixed with
  mhcgnomes' own `species="HLA"` strict constraint on `parse_mhc`, which
  restores the HLA-only contract without reintroducing the lowercase bug —
  confirmed empirically for SLA/BoLA/DLA/RT1 designations and for `a2`.

`parse_mhc` itself is now defensive end to end: a non-string `value`, an
unrecognized `expect`, or any exception raised inside mhcgnomes (including
the deferred import) all degrade to `None` per its documented contract,
rather than raising `AttributeError`/`KeyError`/an uncaught exception. That
also closes the one other real gap the review found: `normalize_mhc_
restriction` no longer has (or needs) its own separate exception guard,
because the function it calls now genuinely never raises.

The remaining findings were cleanup: `_EXPECTED_RESULT_TYPES`'s
`import_module`/`getattr` indirection (the direct cause of the `KeyError`
footgun) is gone in favor of plain deferred class imports; `_parse_hla`'s
dead one-line passthrough is removed and `normalize_mhc_restriction` calls
`parse_mhc` directly; the duplicated `HLA-` prefix-stripping in
`spanning._allele_locus` now calls the newly-public `mhc.strip_hla_prefix`
instead of hand-rolling it a second time; `serotype_key`'s docstring no
longer claims it resolves a split serotype to its broad parent (verified
`serotype_key('A2403') == 'A2403'`, not `'A24'` — that resolution is
`serotype_keys` reading both names hitlist stores in the `serotypes` cell,
not anything `serotype_key` does alone); and `_filter_by_serotype`'s
row-wise `.map(lambda cell: ...)` is now map-over-uniques-then-broadcast,
matching the reviewer's benchmarked ~6x finding.

24 new/extended tests across `test_mhc.py`, `test_cli_hits.py`, and
`test_spanning.py`. Required checks: `./format.sh` (one file reformatted),
`./lint.sh` (clean), `./test.sh` — 478 passed, 1 pre-existing failure
(`test_global53_default_uses_mhcflurry_runtime_calibration_when_available`,
confirmed identical on unmodified main: a shared-environment mhcflurry
version drift unrelated to this diff, not one of the 11 findings).

Version bumped 1.25.4 -> 1.25.5.

## #146 pre-merge review

Only action references and the version changed outside this log. Published
metadata confirms Node.js 24 for checkout v7, setup-python v7, deploy-pages
v5, and upload-pages-artifact v5's nested upload-artifact v7. All workflows
parse, and a semantic before/after comparison proves their non-version
configuration unchanged. Format/lint passed; full suite: 561 passed,
7 warnings. The new workflow will be exercised by the PR CI matrix.

## Release verification: #132

PR #167 merged at 0b28a31; issue #132 is closed. `./deploy.sh` passed
lint and all 561 tests from clean main and published 1.31.1. PyPI reports
both wheel and sdist, with SHA-256 digests matching the local build files.
#168 has passed its full Python 3.9–3.12 CI matrix with no Node deprecation
annotation (only the unrelated upcoming Ubuntu runner-image notice).

## #160 pre-merge review

Eight new tests failed on the original implementation; all ten regression
cases now pass. Existing flagged-target behavior and its rationale/overlap
warning remain intact, while a strict gene cannot evade confidence or mTEC
selection by changing category. mTEC and measured tumor TPM apply to both
categories. Diagnostics distinguish unknown symbols from curated exclusions.
The full suite passed: 571 tests, 7 warnings; format and lint passed.
Clinical documentation correction is tracked in #170; FDA confirms afami-cel
is MAGE-A4-directed, and NY-ESO-1 trial evidence is cited separately.

## #130 pre-merge review

All ten HPA v23 protein-expression/reliability values were recomputed from
source and matched the current oncoref table. The published 1.8.150 minimum
wheel also has those corrected values, so no dependency bump or local
regeneration is appropriate. Ten regression cases now guard these measured
IHC calls and expressed status without changing canonical admission decisions.
Format/lint passed; the full suite passed 581 tests (17 warnings).

## #131 pre-merge review

Current canonical data already satisfy SUN5/SUN3 seed inclusion. Five tests
now guard their identities, RNA/protein distinction, and selection boundary,
including keeping SUN1/SUN2/SPAG4 outside default selection. Re-read the
SUN5 primary paper and human SPAG4 reports; HPA v23 directly verifies SPAG4
pancreas 87.3 and testis 31.1 nTPM. Candidate-reference discoverability is
tracked upstream in oncoref#548. No tissue-safety inference is weakened.
Format/lint passed; full suite: 586 passed, 19 warnings.

## #120 pre-merge review

Documented the current HERV-K/HML-2 coverage limitation at all relevant
entry points. Read Telescope and ERVmap primary methods and HERV-K Env CAR
experimental evidence; the text distinguishes gene/family/locus expression,
translation, presentation, and safety. No nonexistent adapter is advertised.
Strict MkDocs build passed, including the new anchors. Format/lint passed;
full suite after restoring the current reference: 586 passed, 19 warnings,
no skips. Version is 1.31.6.

## Release checkpoint before the final merge (2026-09-22)

| Issue | Reviewed PR | Version | Verified state |
|---|---|---|---|
| #132 | #167 | 1.31.1 | Merged; wheel/sdist published and hashes verified |
| #146 | #168 | 1.31.2 | Merged; wheel/sdist published and hashes verified |
| #160 | #171 | 1.31.3 | Merged; wheel/sdist published and hashes verified |
| #130 | #172 | 1.31.4 | Merged; wheel/sdist published and hashes verified |
| #131 | #173 | 1.31.5 | Merged; wheel/sdist published and hashes verified |
| #120 | #174 | 1.31.6 | Reviewed and authorized; final log update awaits CI/merge/deploy |

During the session, another development update advanced pyensembl to 2.10.17
and selected the full patch/haplotype annotation. Indexed that selected GTF
without deleting the old cache, then reran the complete final suite: 586
passed with no skips. The real gene-length and fragment-model tests passed.
Revalidated local import paths for all 14 relevant packages; changes from
other current development sessions (including hitlist 1.62.21 and varcode
9.4.0) are exercised directly rather than replaced with release wheels.

Upstream issues: openvax/sercol#4 (existing metadata conflict),
openvax/topiary#374 (new osteosarc pin report), pirl-unc/oncoref#548
(SPAG4 candidate-reference omission). Local review findings #169 and #170
are fixed by #171. Existing oncoref#544 remains the separate CTAG2
normal-tissue investigation; no safety conclusion was invented here.

The user directly approved all remaining merges and PyPI deployments. Each
of 1.31.2–1.31.5 was published with `./deploy.sh` from clean main before the
next PR merged. Their deployment suites passed 561, 571, 581, and 586 tests,
respectively. Downloaded wheels and sdists match the corresponding local
builds byte-for-byte and match PyPI's SHA-256 metadata. Each merged PR now
contains its merge commit, test result, release link, and artifact digests.
The final post-merge 1.31.6 audit will be recorded on PR #174 to avoid claiming
publication before it occurs.

Rechecked all 14 local development import paths during the release sequence;
this also exercised current topiary 5.68.2 and varcode 9.4.1. A clean archive
build of the 1.31.6 source passed the version/source/package-data audit for
both wheel and sdist. No duplicate CTA definition/proteoform tables or
regeneration sidecars are packaged.

### Next work, ordered by dependency and urgency

1. Resolve development-environment metadata conflicts in
   [sercol#4](https://github.com/openvax/sercol/issues/4) and
   [topiary#374](https://github.com/openvax/topiary/issues/374); also audit
   [hitlist#504](https://github.com/pirl-unc/hitlist/issues/504)'s required
   import behind an optional dependency declaration. These unblock reliable
   installation and integration across multiple consumers.
2. Settle the full-transcriptome versus gene-subset clean-TPM contract in
   [oncoref#541](https://github.com/pirl-unc/oncoref/issues/541). Its updated
   investigation identifies the pan-cancer subset input as the problem;
   do not change valid full-transcriptome normalization based on its old title.
3. Continue evidence interpretation in
   [oncoref#544](https://github.com/pirl-unc/oncoref/issues/544) (CTAG2 cardiac
   expression) and [hitlist#490](https://github.com/pirl-unc/hitlist/issues/490)
   (queried versus observed allele scope); retain SPAG4 candidate-reference
   follow-up in [oncoref#548](https://github.com/pirl-unc/oncoref/issues/548).

## MHCflurry panel contract audit (2026-09-29)

### Specification

Establish the exact software version, model download release, named sequence
columns, full pseudosequences, and calibration behavior behind the C*14:02 /
C*14:03 claim. Separate representation changes from model training evidence and
from tsarina's stable, named panel membership. Fix the demonstrated integration
failure without inventing biological equivalence or changing a 53-allele panel
into 54 alleles under its existing name. Audit proposed dependency-floor and
pandas changes against actual requirements rather than release age.

- [x] Reproduce the failing calibration test and compare old/new model artifacts.
- [x] Inspect training provenance for C*14:03 and record reproducible evidence.
- [x] File confirmed problems on their owning repositories; link from the PR.
- [x] Implement the justified panel/test/CI contract and document model scope.
- [x] Repair local editable installs with `./develop.sh` and verify import paths.
- [x] Run `./format.sh`, `./lint.sh`, and `./test.sh`; review the complete diff.
- [ ] Bump version on the feature branch, open/merge PR, deploy clean main.
- [ ] Verify PyPI artifacts and review the next dependency-ordered issue block.

### Plan check-in

The implementation follows the artifact audit. Dependency floors will change
only if a required API/data contract demonstrably needs the newer minimum;
a pandas upper bound requires an observed incompatibility. Named panel changes
need evidence and an explicit compatibility policy.

### Review

PR: https://github.com/pirl-unc/tsarina/pull/183 (fixes #181 and #182).

- `./format.sh` and `./lint.sh`: passed.
- `./test.sh --run-mhcflurry`: 609 passed, no skips, pandas 2.3.3 and
  current editable sibling dependencies; 81% coverage.
- Same full suite with isolated pandas 3.0.6 and published oncoref 1.8.206:
  609 passed, no skips; 81% coverage. No pandas upper bound is warranted by
  this suite. This is tsarina compatibility evidence, not exhaustive upstream
  pandas validation.
- Real-model integration against the legacy 2.2.0 bundle: 2 passed.
- Missing-model negative check: enabled integration fails explicitly.
- QC negative check: injecting GAGE12B into canonical targets fails explicitly.
- Strict MkDocs build passed; CLI default/legacy override verified.
- CI at 0266798 passed lint, Python 3.9–3.12 unit tests and real-model integration.

The shared environment loads unrelated third-party pytest plugins. Final local
runs used `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1` with explicit pytest-cov/xdist
plugins after interrupting slow concurrent runs; no tests were excluded.
Deployment uses the same explicit plugin configuration and enables model tests.
Final merge/publication evidence will be recorded on PR #183 after execution,
avoiding a commit on main or a claim of publication before upload succeeds.

### Dependency and QC findings

The published oncoref 1.8.206 wheel and the current checkout both include
rejected/noncoding candidates in raw evidence. All nine protein-model flags
are outside `cta_gene_ids()` (published wheel: 2,532 evidence rows / 624
canonical IDs). Issue #182 tracks the incorrect raw-universe assertion.
The test now checks canonical admitted IDs while retaining the complete raw
diagnostic. No membership, evidence, or length threshold changes.

The pandas 3.0.6 baseline passed 603 tests and failed only that same QC
assertion; it also failed on pandas 2.3.3. No pandas-specific incompatibility
supports the proposed upper bound. Current floors retain their audited API
contracts rather than being raised to whichever release happens to be installed.

Real-model calibration and scoring passed against the 2.3.0 presentation bundle.
All seven old panel memberships compare exactly equal to origin/main; only the
new global54 panel and CLI/API defaults change. The release version is 1.32.0.
