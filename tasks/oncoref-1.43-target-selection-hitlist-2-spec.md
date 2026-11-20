# OncoRef 1.43.0 Rebuilt Target Selection & Hitlist 2.1.22 EvidenceStore Adoption Spec

Date: 2026-11-20
Branch: `adopt-oncoref-1.43-hitlist-2.1.22`
Target Tsarina version: `1.39.0`

## 1. Goals

1. **Adopt OncoRef 1.43.0's rebuilt set-based target catalog (`oncoref.target_selection` / `oncoref.antigen_sets`)** as Tsarina's default human CTA definition (`oncoref.target-selection.2026-11-19.4`).
2. **Adopt Hitlist 2.1.22** as Tsarina's MS evidence layer:
   - Integrate Hitlist 2.0's verified `EvidenceStore` packages as a first-class runtime MS backend alongside `observations.parquet`.
   - Regenerate `tsarina/data/gene-ms-safety-evidence.csv` and `tsarina/data/gene-ms-safety-evidence.meta.json` from a Hitlist 2.1.22 corpus covering every current `target_selection` CTA, candidate (`CA|RA`), and previously tracked gene.
   - Align Python and dependency floors (`requires-python >= 3.11`, `oncoref >= 1.43.0`, `hitlist >= 2.1.22`).

## 2. OncoRef 1.43.0 `target_selection` Contract

- **Policy**: `oncoref.target-selection.2026-11-19.4` (19,893 protein groups / 20,099 member genes).
- **Categories**:
  - `CTA`: 327 protein groups / 358 member genes (`is_cta == True`), with `cta_tier` in `("core", "core_with_extended_exposure", "extended_only")` (206 / 100 / 21).
  - `Cancer Antigen only`: 70 protein groups (`is_cancer_antigen == True`, `is_cta == False`).
  - `Reproductive Antigen only`: 128 protein groups (`is_reproductive_antigen == True`, `is_cta == False`).
  - `Other`: 1,118 protein groups (`Other:Review` 734, `Other:Expression` 248, `Other:Unscreened tissue` 88, `Other:Sequence` 48).
  - `Unsafe`: 18,250 protein groups.
- **Scientific status disclosure**: `target_selection.metadata()["release_ready"]` is `False` while three open research reviews remain (tumor-source comparability, read-level rare-cell mapping-risk confirmation, hair-follicle normal-protein replay #687); unresolved evidence stays excluded from positive sets. Tsarina propagates `cta_selection_release_ready` and `cta_selection_open_reviews` in `selection_provenance()`.
- **Gene identity & compatibility**:
  - Canonical symbols, synonyms, and identity status come from `oncoref.gene_identity_catalog()` (20,097/20,099 members; the 2 non-catalog loci keep `gene_id` as symbol and `identity_status="registry_member_outside_canonical_catalog"`).
  - `CTA_selection(category=None, *, species="human")` returns a detached normalized DataFrame from `oncoref.target_selection.results()` for human (with `gene_ids`, `gene_identities`, `policy_id`, `reference_id`, `source_snapshot_id`, `reasons`, `evidence_flags`, `cta_tier`, `classification`, `is_cta`, `is_cancer_antigen`, `is_reproductive_antigen`, `is_target_candidate`), and delegates non-human species to `oncoref.cta_selection.results()`.
  - `CTA_evidence()` retains `selection_*` prefixed columns (including `selection_cta_tier`, `selection_classification`, `selection_is_cancer_antigen`, `selection_is_reproductive_antigen`, `selection_is_target_candidate`), never adding unprefixed `cta_tier`.
  - `CTA_gene_ids/names()` = `CTA` (358 genes); `CTA_excluded_gene_ids/names()` = `Unsafe`; `Cancer Antigen only`, `Reproductive Antigen only`, and `Other` remain distinct and stay in the non-CTA sequence background.
  - New explicit helpers: `CTA_by_tier(tier)`, `CTA_cancer_antigen_gene_ids/names()`, `CTA_reproductive_antigen_gene_ids/names()`, `CTA_target_candidate_gene_ids/names()`.
  - `load_vaccine_inputs(definition="current")` consumes the rebuilt `target_selection` catalog and records its policy, reference, snapshot, and `release_ready` receipts.

## 3. Hitlist 2.1.22 MS Data & `EvidenceStore` Contract

- **New module `tsarina/evidence_store.py`**:
  - `resolve_evidence_package(package_path=None)`: returns explicit path or `TSARINA_EVIDENCE_PACKAGE` / `HITLIST_EVIDENCE_PACKAGE` env var.
  - `open_evidence_store(package_path)`: cryptographically verifies and caches `hitlist.EvidenceStore.open(path)` keyed by `(resolved_path, manifest_mtime_ns, manifest_size)`.
  - `evidence_store_provenance(package_path)`: returns manifest SHA256, package/catalog schema versions, Hitlist versions, and `mhc_positive_ms` count.
  - `load_evidence_store_ms_hits(package_path, *, peptides=None, mhc_class="I", mhc_species="Homo sapiens", gene_name=None, columns=None, include_binding=False, include_other=False, exclude_non_peptide_ligand=True, length_min=None, length_max=None)`:
    - Queries `mhc_positive_ms` (plus optional `mhc_binding` / `mhc_other`) through `EvidenceStore.query(...)` with Arrow pushdown filters.
    - Resolves `gene_name` via the verified package's `native/peptide_mappings.parquet` when needed.
    - Projects deferred native fields (`mhc_class`, `mhc_species`, `species`, `allele_resolution`, `serotypes`, `serotype_source`, `source_tissue`, `cell_name`, `cell_line_name`, `src_*`, `is_monoallelic`, `mhc_allele_provenance`, `mhc_allele_set`) directly from the manifest-verified `native/<partition>` by `native_row_index`, without scanning the entire native table.
    - Returns a DataFrame compatible with `load_observations` consumers plus typed fields (`observation_id`, `outcome_status`, `restriction_evidence`) and attaches cryptographic package provenance to `df.attrs`.
- **Wiring**:
  - `tsarina.indexing.load_ms_evidence(..., evidence_package=None)` routes through `load_evidence_store_ms_hits` when `evidence_package` (or env var) is configured; otherwise queries Hitlist 2.1.22 `load_ms_observations` with `exclude_non_peptide_ligand=True`.
  - Thread `evidence_package` through `load_public_ms_hits`, `cta_healthy_tissue_ms_hits`, `CTA_detailed_evidence`, `personalized_targets`, `target_peptides`, `cli_hits` / `cli_common` (`--evidence-package`).
- **Bundled overlay regeneration**:
  - `scripts/regenerate_gene_ms_safety_evidence.py` searches all current `CTA` + `CA|RA` candidate genes plus previously tracked genes, supports both `--evidence-package` and default `observations.parquet`, and writes `gene-ms-safety-evidence.csv` + `gene-ms-safety-evidence.meta.json`.
  - Rebuild an isolated Hitlist 2.1.22 observations index under `/tmp/tsarina-hitlist-2.1.22` (never modifying `~/.hitlist`) from the local IEDB + CEDAR exports and regenerate the bundled overlay.

## 4. Verification Criteria

- Every `target_selection` CTA gene (358) and `CA|RA` candidate gene (562) is searched in `gene-ms-safety-evidence.csv` and verified against `gene-ms-safety-evidence.meta.json`.
- Real `hitlist.write_mhc_evidence_package` -> `EvidenceStore.open` -> `load_evidence_store_ms_hits` / `load_ms_evidence` / `load_public_ms_hits` / `cta_healthy_tissue_ms_hits` / `tsarina hits --evidence-package` end-to-end tests pass (including tampered-manifest rejection!).
- `./format.sh`, `./lint.sh`, `./test.sh --run-mhcflurry`, `mkdocs build --strict`, and `python -m build && twine check dist/*` all pass.
- Clean PR, CI, merge, and PyPI deploy (if network is available; if offline, document exact state and next step clearly).
