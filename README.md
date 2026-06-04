# Kpneumoniae-ST11-Plasmid-Network
# A dual‑tier plasmid network drives the evolutionary success of a pandemic *Klebsiella pneumoniae* lineage

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Manuscript under review at *Scientific Reports*.** A DOI will be linked here upon availability.

---

## Overview

This repository contains the complete computational workflow used for our study of the genomic epidemiology and plasmid network architecture of *Klebsiella pneumoniae* ST11.  
We integrated pan‑genomics, phylogenetics, plasmid co‑occurrence network analysis, and statistical enrichment testing to dissect the evolutionary mechanisms underlying the pandemic success of the ST11 lineage.

**Key findings include:**

- **Clonal sweep** – ST11 accounts for 30% of the global collection, characterized by a severe phylogenetic bottleneck and a ~20-fold reduction in Mean Phylogenetic Distance (MPD = 0.0049) compared to the species-wide background.
- **Pangenome dynamics** – The species‑wide pangenome is extremely open (Heaps' law α = 0.59), while the East Asian ST11 clade exhibits the characteristically low diversity expected of a recent clonal expansion (α = 0.87).
- **Plasmid network architecture** – Lineage‑specific “plasmotypes” and hub replicons (IncFII, IncFIB) form a conceptual dual‑tier network model that balances vertical stability and horizontal gene transfer.
- **Chromosomal integration** – Tandem amplification and chromosomal capture of key resistance determinants (*bla*<sub>KPC‑2</sub>, *bla*<sub>CTX‑M‑15</sub>) reinforce stable foundational resistance.
- **Convergence** – 34.2% of ST11 isolates carry both carbapenem resistance and hypervirulence markers.

---

Most scripts are self‑contained; intermediate outputs are written to a `results/` directory (not included in this repository due to size).

---

## Requirements

### Command-line tools

- **NCBI Datasets CLI** – for downloading genome assemblies from NCBI.
- **jq** – JSON processor (used to parse metadata).
- **GNU Parallel** – for parallel execution of shell commands.
- **mlst** (Seemann T) – multi-locus sequence typing for K. pneumoniae.
- **Kleborate v3.2.4** – genome‑based ST, K‑type, and O‑type assignment; includes Kaptive for capsule typing.
- **Panaroo v1.5.2** – pangenome construction and gene presence/absence analysis.
- **IQ‑TREE v3.0.1** – maximum‑likelihood phylogenetic tree inference.
- **iTOL v6** – interactive visualization and annotation of phylogenetic trees (used via the web interface).
- **ResFinder v4.7.0** – identification of acquired antimicrobial resistance genes (uses KMA for read mapping).
- **KMA v1.6.8** – k‑mer alignment for mapping short reads or contigs to a reference database.
- **PlasmidFinder v2.1.6** – detection of plasmid replicons (Enterobacteriaceae database from CGE).
- **MOB‑suite v3.1.9** – prediction of plasmid mobility and conjugation potential.

### R packages

- micropan
- picante
- ggplot2
- dplyr
- ape
- data.table
- tidyr

### Python packages

- pandas
- networkx
- matplotlib
- numpy
- scipy
- statsmodels
- seaborn
