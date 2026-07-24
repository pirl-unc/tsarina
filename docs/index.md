# Tsarina documentation

Tsarina turns shared cancer targets into ranked peptide-MHC candidates. Use it
to prioritize CTA, viral, and recurrent-mutation targets for one patient, or to
design an off-the-shelf CTA panel across a population HLA set.

## Start with your outcome

### Prioritize targets for one patient

Use patient HLA type plus any available CTA expression, hotspot mutations, and
viral status. Tsarina returns ranked pMHCs with source abundance, public MS
support, healthy-tissue safety flags, and predicted presentation.

Continue to [Personalized target selection](personalized-targets.md).

### Design a reusable CTA panel

Use Tsarina's automatic CTA safety filters and cancer-prevalence ranking, or
provide an explicit CTA list. Tsarina returns a CTA × HLA matrix, evidence
tiers, and population-coverage estimates.

Continue to [CTA panel design](panel-design.md).

### Inspect peptide evidence

Use the data registry to install or register IEDB and CEDAR exports, then query
the observation index for specified peptides.

Continue to [Data and evidence](data-and-evidence.md).

## The shared pipeline

All target-selection workflows follow the same conceptual stages:

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

These stages keep target definition distinct from downstream evidence. In
particular, [oncoref](https://github.com/pirl-unc/oncoref) alone owns CTA
membership, aliases, HPA restriction calls, and proteoform groups. See
[CTA ownership and downstream evidence](curation.md) for the integration
boundary.

## Install and prepare data

Install the package:

```bash
pip install "tsarina[all]"
```

See which external datasets Tsarina can use:

```bash
tsarina data available
```

Viral proteomes can be fetched automatically. IEDB and CEDAR ligand exports
must be downloaded under their source terms and registered locally:

```bash
tsarina data fetch hpv16
tsarina data register iedb /data/mhc_ligand_full.csv
tsarina data register cedar /data/cedar-mhc-ligand-full.csv
tsarina data list
```

The [Data and evidence](data-and-evidence.md) guide explains storage,
observation classification, and the evidence model.

## Core target concepts

### Cancer-testis antigens

CTAs are normally restricted to reproductive tissues but can be reactivated in
tumors. Oncoref supplies the candidate universe and HPA-derived restriction
axes. Tsarina adds tumor prevalence, public-MS safety evidence, peptide
enumeration, and HLA scoring after membership is known.

### Oncogenic viruses

Viral proteins are foreign to the human proteome. Clinical helpers default to
human-exclusive viral peptides, removing viral k-mers that also occur in human
proteins.

### Recurrent mutations

Hotspot driver mutations produce shared mutant peptides in patients carrying
the same alteration. Tsarina enumerates mutation-spanning peptides rather than
discovering private passenger mutations from whole-exome sequencing.

## Guide map

| Guide | Use it for |
|---|---|
| [Personalized target selection](personalized-targets.md) | Patient inputs, CLI and Python use, output columns, and ranking |
| [CTA panel design](panel-design.md) | Automatic CTA selection, safety gates, evidence tiers, HLA panels, and coverage |
| [Data and evidence](data-and-evidence.md) | Dataset setup, source classification, positive/negative evidence, scoring, and naming |
| [CTA ownership and downstream evidence](curation.md) | The oncoref/Tsarina boundary, API behavior, and definition updates |

For the current command-line surface, use:

```bash
tsarina --help
tsarina personalize --help
tsarina panel --help
tsarina hits --help
tsarina data --help
```
