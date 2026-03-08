# Saccharomyces Phylogenetic Network — Project Brief

## Goal

Build a phylogenetic network (reticulate DAG) of the *Saccharomyces* genus, covering all eight *sensu stricto* species, their known interspecific hybrids, and the major strain-level populations of *S. cerevisiae* relevant to fermentation. No such database currently exists in graph form; the raw scientific data is available but has never been consolidated into a queryable network structure.

## Representation decisions

- **Primary structure**: a Python `NetworkX` directed graph (`DiGraph`), with hybrid/reticulation nodes having two incoming edges (one per parent).
- **Serialisation format**: Extended Newick (Rich Newick) for interoperability with phylogenetics tools (PhyloNet, IcyTree, Dendroscope, SplitsTree). The `#H` label convention marks hybrid nodes.
- **Two-level model**:
  - **Level 1 (species network)**: ~10 nodes, ~5 hybrid nodes. Fully tractable, well-supported by literature. Build this first.
  - **Level 2 (strain/population network)**: extend Level 1 with *S. cerevisiae* population nodes (Beer 1, Beer 2, Wine/European, Sake, West African, Malaysian, North American, farmhouse, wild) and named hybrid strains. More nodes, same graph structure.

## Node schema

Each node should carry at minimum:

```python
{
    "id":            str,       # e.g. "S_cerevisiae", "S_pastorianus"
    "type":          str,       # "species" | "hybrid" | "population" | "strain"
    "display_name":  str,
    "is_hybrid":     bool,
    "fermentation":  list[str], # e.g. ["ale", "wine", "sake"] or [] for wild
    "wild":          bool,
    "geography":     list[str], # known natural ranges
    "notes":         str
}
```

## Edge schema

```python
{
    "type":               str,    # "divergence" | "hybridisation" | "introgression"
    "admixture_fraction": float,  # 0.0–1.0; 0.5 for equal hybrid, None for pure divergence
    "divergence_mya":     float,  # millions of years ago, where known; None otherwise
    "confidence":         str     # "high" | "medium" | "low"
}
```

## Level 1 species tree topology

The backbone tree (before adding reticulations) is:

```
Outgroup (Naumovozyma castellii)
└── Saccharomyces sensu stricto
    ├── S. arboricola          (East Asian oaks; wild only)
    └──┬── S. kudriavzevii     (European/Asian oak bark; wild)
       └──┬── S. mikatae       (Japan; wild only)
          ├── S. jurei         (European pre-Alps/UK oaks; wild; described 2017)
          └──┬── S. paradoxus  (global; wild, rare spontaneous fermentation)
             └──┬── S. uvarum  (wine, cider, cold fermentation; partly domesticated)
                └──┬── S. eubayanus  (Patagonia/Tibet; wild; lager parent)
                   └── S. cerevisiae (global; primary domesticated fermentation yeast)
```

Divergence timescale: the entire sensu stricto crown group is ~20 million years old. *S. cerevisiae* / *S. eubayanus* split is the most recent speciation among the eight species.

## Known interspecific hybrid nodes (reticulations to add)

| Hybrid node | Parent 1 | Parent 2 (and 3) | Context | Admixture | Confidence |
|---|---|---|---|---|---|
| *S. pastorianus* | *S. cerevisiae* | *S. eubayanus* | Lager beer; two independent hybridisation events (Saaz/Group 1 and Frohberg/Group 2) | ~0.5 each, but Group 1 is ~haploid for Sc, Group 2 ~diploid | High |
| *S. cerevisiae* × *S. kudriavzevii* hybrids | *S. cerevisiae* (wine strains) | *S. kudriavzevii* (European) | Central European wine fermentation; at least 6 independent events inferred | Variable; proportion of Sk genome differs per strain | High |
| *S. cerevisiae* × *S. uvarum* hybrids | *S. cerevisiae* | *S. uvarum* | Wine and cider | Variable | Medium |
| *S. cerevisiae* × *S. kudriavzevii* × *S. uvarum* triple hybrids | *S. cerevisiae* | *S. kudriavzevii* + *S. uvarum* | Wine | Variable | Medium |
| *S. bayanus* (CBS 380T) | *S. uvarum* (~67%) | *S. eubayanus* (~33%) + *S. cerevisiae* introgression | Historical type strain; effectively a three-way mosaic | ~0.67 / ~0.33 | High |

Note: *S. kudriavzevii* has never itself been isolated from fermentative environments; the hybridisation is believed to have occurred in the wild before hybrid strains colonised fermentation niches.

## Level 2: S. cerevisiae population nodes

These sit below the *S. cerevisiae* species node as a subgraph. Populations from Gallone et al. (2016) and Liti et al. (2009):

| Population node | Fermentation association | Geography | Notes |
|---|---|---|---|
| Beer 1 | Ale (British, Belgian abbey, many craft strains) | Europe/global | Largest domesticated ale cluster |
| Beer 2 | Belgian saison, wheat beer | Belgium/Europe | Distinct from Beer 1 |
| Wine / European | Wine, bread | Europe, global | Largest overall; includes many commercial wine yeasts |
| Sake | Sake | Japan | |
| West African | Palm wine, local fermentation | West Africa | |
| Malaysian | Wild/fermentation | Southeast Asia | |
| North American | Wild | North America | |
| Farmhouse (European landrace) | Farmhouse ale (kveik, Lithuanian, Latvian) | Norway, Baltics | Mixed Beer 1 + Asian domesticated ancestry; Preiss et al. 2024 |
| Wild (root lineages) | None | Global | Oldest lineages; highest genetic diversity; at root of S. cerevisiae tree |
| Andean chicha | Chicha (maize beer) | Ecuador, Peru | Related to Mexican agave and French Guiana strains; carries STA1 diastatic gene |

*S. cerevisiae* also exhibits within-species introgression from *S. paradoxus* (notably enriched in Neotropical strains — Mexico, French Guiana, Ecuador, Brazil). This can be represented as a low-weight introgression edge from *S. paradoxus* to specific population nodes.

## Primary data sources

| Source | What it provides | Access |
|---|---|---|
| Peris et al. (2023), *Nature Communications* 14:690 | Geographies, hosts, substrates, phylogenetic relationships for ~1,800 strains; 163 complete genomes; 128 phenotyped strains | https://www.nature.com/articles/s41467-023-36139-2 |
| Sac2.0 companion site | Genome assemblies (ENA: PRJEB48264), raw reads (SRA: PRJNA475869), BLAST server, interactive tree | https://perisd.github.io/Sac2.0 / https://github.com/PerisD/Sac2.0 |
| Gallone et al. (2016), *Cell* 166:1397 | Ale yeast Beer 1 / Beer 2 classification; strain-level dendrogram of domesticated S. cerevisiae | https://doi.org/10.1016/j.cell.2016.08.020 |
| Liti et al. (2009), *Nature* 458:337 | Five S. cerevisiae populations (North American, Malaysian, West African, Sake, Wine/European); mosaic strains | https://doi.org/10.1038/nature07743 |
| Milk the Funk wiki | Aggregated summary of brewing strain phylogenetics; maps Gallone codes to commercial White Labs/Wyeast products | https://www.milkthefunk.com/wiki/Saccharomyces |
| ScRAPdb (2025), *Nucleic Acids Research* | Pan-omics database for S. cerevisiae strains; interactive phylogeny, pangenome | https://doi.org/10.1093/nar/gkae1020 |

## Suggested build order

1. Implement the Level 1 species network in NetworkX (8 species + outgroup + 5 hybrid nodes).
2. Validate topology against Extended Newick export; visualise with IcyTree or Dendroscope.
3. Extend to Level 2 by adding *S. cerevisiae* population nodes as children of the species node.
4. Add introgression edges (e.g. *S. paradoxus* → Neotropical *S. cerevisiae* populations).
5. Populate node metadata from Peris et al. supplementary data (strain counts, geographies, fermentation associations).
6. Export to both Extended Newick and GraphML for downstream use.

## Key open question deferred

Whether to extend to a full Ancestral Recombination Graph (ARG), where admixture proportions vary per genomic locus, was discussed and set aside for now. The strain-level data in Peris et al. and ADMIXTURE analyses from Gallone et al. would be the input if this is revisited. Tools: `tsinfer`, `Relate`, `msprime`.
