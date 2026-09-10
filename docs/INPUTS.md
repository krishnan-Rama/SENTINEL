# Inputs

SENTINEL needs three things from you. Everything else it fetches or computes.

## 1. A cleaned toxicity table

One row per species, tab-separated, pointed at by `ECOTOX_TABLE` in the config. The
reference dataset was built from the US EPA ECOTOX knowledgebase.

| Column | Meaning |
|--------|---------|
| `Species_Name` | binomial, as it should appear in outputs |
| `Class` | taxonomic class, used for the clade decomposition |
| `Conc_Mean_Std` | the endpoint on a mass basis, mg/L |
| `Conc_Type` | active ingredient or formulation |
| `Pub_Year` | earliest publication year contributing to the value |

Curation decisions to make before this file exists, because the pipeline cannot make them
for you: which exposure duration and endpoint type to accept, whether to pool formulations
with technical-grade material, how to combine multiple records for one species, and what to
do with records that resolve only to genus or family. Record what you chose. The reference
dataset used 96 h mortality LC50 values on an active-ingredient basis from 1990 onwards,
geometric means where a species had several records, and dropped records resolving only to
family.

The pipeline tests whether the species you lose to missing sequence data differ in
sensitivity from the ones you keep. That test is only meaningful if this table is the full
curated set, before any sequence availability was considered. Do not pre-filter it.

## 2. A contact address

`UNIPROT_CONTACT`. Sent in the User-Agent header of UniProt requests. Requests without one
are throttled harder and the maintainers ask for it.

## 3. A dated species tree, if you want the analysis to be publishable

`SPECIES_TREE`, a Newick file whose tip labels match the species names in your endpoint
table. TimeTree exports work after pruning to your species set.

Without it the pipeline falls back to the target gene tree. That is adequate for deciding
what a sequence is and inadequate for anything phylogenetic downstream, for two reasons: a
gene tree of a duplicated family groups by paralogue before taxonomy, and branch lengths
measure sequence divergence rather than time. On the reference dataset the two choices gave
Pagel's lambda of 0.885 and 0.965 on the same endpoints, which is the difference between
"strong signal" and "almost entirely phylogenetic".
