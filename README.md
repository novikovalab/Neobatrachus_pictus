Analysis scripts for the study:
“Changes in cell cycle machinery as genetic basis for polyploidy stabilization in Australian burrowing frogs (Neobatrachus)”

This repository contains main computational and bioinformatics scripts used in the study
“Changes in cell cycle machinery as genetic basis for polyploidy stabilization in Australian burrowing frogs (Neobatrachus)”.
The scripts span the full analysis workflow used in the manuscript, including:

population genomic analyses (PCA, admixture, NJ trees, SNP sharing),

mixed-ploidy network and inheritance mode inference,

gene family evolution (CAFE5),

genome synteny and dotplot visualization,

heatmap-based selection scan summaries,

stLFR/Supernova assembly pipelines,

and protein structural modelling.

The scripts are provided as standalone files and can be adapted to similar datasets or species.

📂 Contents
Population genomics & evolutionary analyses

PCA admixture tree barplot.R – PCA, distance matrices, hierarchical clustering, NJ trees, and admixture plotting.

Run_network_PCA_tree.R – PCA-based mixed-ploidy network reconstruction.

Run_shared_uniq_SNPs.R – Identification and visualization of shared and lineage-specific SNPs across genomes.

Run_inheritance_mode.R – Haplotype assignment and inheritance mode analysis for tetraploid species.

cafe_analysis – CAFE5 analyses for gene family expansion / contraction across species.

Genome assembly & synteny

Run_CoGe_1.sh – First-step processing of CoGe synteny outputs.

Run_CoGe_2_synteny block processing + circular plot.R – Synteny block parsing and circular ideogram visualization.

Run_assembly_N. kunapalari, N. sudellae, and others – stLFR Supernova genome assembly and scaffolding pipelines.

Selection scans

selection_scan_genes_heatmap.py – Multi-page heatmap visualization combining multiple gene-level scan results; includes hierarchical clustering and species ordering.

Protein structure modelling

protein modelling pymol.r – PyMOL automation script for modelling the SYCE2–TEX12 complex.

⚠️ Important notes
🔧 1. Some scripts are merged from multiple original scripts

Several scripts in this repository combine multiple analysis steps into a single file (e.g., synteny processing, PCA + admixture + tree plotting, or multi-gene heatmaps).
These integrated versions simplify the workflow but may differ from the original pipeline structure described in internal notes.

📁 2. File paths are not included and must be adapted

For security and portability, explicit file paths (e.g., /netscratch/...) have been removed or replaced with placeholders.
Users must modify:

input and output directories

reference genome locations

VCF or matrix file paths

synteny output files

annotation / GFF / BED resources

based on their own environment and data organization.

🧬 3. Species lists must be adjusted manually

Several scripts rely on species order, ploidy assignment, or naming patterns (e.g., PCA coloring, network analyses, SNP-sharing functions).
Users working with different species, sample sizes, or naming conventions must modify:

species vectors

color palettes

ploidy mappings

ordering of samples in PCA/NJ/admixture plots

patterns used for grouping samples

🌱 4. External software and environment dependencies are not bundled

The repository does not include the full list of software modules or exact versions used on the HPC cluster.
Users must ensure that required dependencies are available, including:

R packages:
ggtree, factoextra, ape, karyoploteR, GenomicRanges, dendextend, RColorBrewer, etc.

Python packages:
numpy, pandas, matplotlib, scipy, vcfpy or vcfR replacement as needed.

External tools:
CAFE5, CoGe, MAFFT, IQ-TREE, stLFR Supernova, PyMOL, samtools, bcftools, htslib, etc.

Because the HPC environment differs across institutes, users should adapt module loading and paths accordingly.
