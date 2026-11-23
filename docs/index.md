# Tsarina documentation

Tsarina ranks peptide-MHC candidates from cancer-testis antigens (CTAs),
curated cancer fusions, oncogenic viruses, and recurrent hotspot mutations. It
supports single-patient target selection (human and mouse), population CTA
panels, public MS safety lookup, and shared CTA vaccine constructs.

## Workflows

### [Ecosystem and Python API](ecosystem-and-api.md)

How Tsarina connects to `oncoref`, `hitlist`, `pyensembl`, `orfdata`,
`mhcgnomes`, and MHC predictors (`mhcflurry`, `mhctools`, `topiary`), plus a
grouped Python API and CLI command map.

### [Personalized target selection](personalized-targets.md)

Patient HLA (or mouse H-2) type, CTA expression, biological sex eligibility,
hotspot mutations, fusion breakpoints, and viral status determine the
candidates. Results include ranked pMHCs with source abundance, public MS
support, somatic and sex-specific reproductive safety flags, Cys
MS-detectability caveats, and predicted presentation.

### [CTA panel design](panel-design.md)

Automatic selection applies CTA safety filters and cancer-prevalence ranking;
an explicit CTA list is also supported. Results include a CTA × HLA matrix,
evidence tiers and population-coverage estimates.

### [CTA vaccine design](vaccine-design.md)

Rank proteoforms by mortality-weighted p95 prevalence, retain CTA-specific
MS-supported sequence, and audit a single DNA/RNA construct.

See: [Vaccine Atlas](vaccine-results/index.html) — strict/loose antigen
comparisons, source-resolved tissue maps, HLA coverage and cumulative MS evidence.

### [Data and evidence](data-and-evidence.md)

Register IEDB and CEDAR exports, query peptide observations, and interpret
their tissue context, HLA restrictions and prediction scores.

## Target selection

Target selection follows five stages:

1. **Define candidates.** CTA definitions come from oncoref; viral proteins and
   recurrent mutation hotspots come from Tsarina's target modules.
2. **Apply biological context.** Patient workflows retain expressed CTAs and
   detected mutations or viruses. Panel workflows rank CTAs by population-level
   cancer prevalence.
3. **Enforce specificity and safety.** Candidate peptides must be exclusive to
   their intended source rules. Healthy-tissue MS observations are kept
   separate from expected reproductive-tissue and thymus observations.
4. **Add presentation evidence.** Public immunopeptidomics observations and HLA
   presentation predictions establish evidence tiers.
5. **Rank with provenance.** Results retain target identity, observation
   evidence, prediction scores, and the reason for their rank.

[OncoRef](https://github.com/pirl-unc/oncoref) owns CTA membership, aliases,
HPA restriction calls and proteoform groups. Tsarina adds downstream evidence
for selection and ranking.

See: [CTA ownership and downstream evidence](curation.md).

## Install and prepare data

Install the package:

```bash
pip install "tsarina[all]"
```

Inspect local environment readiness and built-in catalogs (read-only):

```bash
tsarina status
tsarina catalog
tsarina data available --grouped
```

IEDB, CEDAR, and built-in viral proteomes can be fetched directly (or
registered locally if already downloaded):

```bash
tsarina data fetch iedb
tsarina data fetch cedar
tsarina data fetch hpv16
tsarina data list
```

See: [Data and evidence](data-and-evidence.md) — storage, observation
classification, preflight checks, and the evidence model.

## Target types

### Cancer-testis antigens (human and mouse)

CTAs are normally restricted to reproductive tissues but reactivated in tumors.
OncoRef supplies current protein-level selection decisions across human (327
CTAs) and evaluated mouse (20 CTAs), biological sex-restricted eligibility
(`--sex female|male`), and historical HPA-derived restriction axes. Tsarina
adds tumor prevalence, public-MS safety evidence (including sex-specific
reproductive flags), peptide enumeration, and MHC scoring.

### Curated cancer fusions

OncoRef curates 188 cancer fusions (`tsarina catalog fusions` /
`cancer_fusions()`). Given explicit fusion breakpoint sequences, Tsarina
enumerates host-exclusive junction-spanning peptides (`fusion_junction_peptides()`
or `--fusion` / `--fusions` in `tsarina personalize`).

### Oncogenic viruses

Viral proteins are foreign to the host proteome. Tsarina supports 9 built-in
UniProt oncogenic viruses plus optional `orfdata` IARC GDC viral subsets
(`cancer-associated`, `cancer-causing`, `cancer-integrated-driver`) and spliced
RefSeq viral ORFs (`tsarina catalog viruses --source orfdata`). Clinical
helpers default to host-exclusive viral peptides.

### Recurrent mutations

Hotspot driver mutations produce shared mutant peptides in patients carrying
the same alteration. Tsarina enumerates host-exclusive mutation-spanning
peptides rather than discovering private passenger mutations from whole-exome
sequencing.

## Guides

| Guide | Use it for |
|---|---|
| [Ecosystem & Python API](ecosystem-and-api.md) | Stack ownership (`oncoref`, `hitlist`, `pyensembl`, `orfdata`, `mhcgnomes`, predictors), CLI map, and Python quickstart |
| [Personalized target selection](personalized-targets.md) | Patient inputs, sex eligibility, fusions, mouse CTAs, output columns, and ranking |
| [CTA panel design](panel-design.md) | Automatic CTA selection, safety gates, evidence tiers, HLA panels, and coverage |
| [Data and evidence](data-and-evidence.md) | Dataset setup, source classification, Cys MS-detectability caveats, scoring, and naming |
| [CTA ownership and downstream evidence](curation.md) | The oncoref/Tsarina boundary, API behavior, and definition updates |

For the current command-line surface, use:

```bash
tsarina --help
tsarina status
tsarina catalog --help
tsarina personalize --help
tsarina panel --help
tsarina hits --help
tsarina data --help
```
