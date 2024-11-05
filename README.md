# SingleCell_Analysis_Scripts

This repository contains R/Python scripts for the analysis of scRNAseq and spatial transcriptomics


## scRNA-seq analysis
Inside the scRNAseq folder, there are scripts covering the following analysis steps:
1. **Quality control** (removal of low quality cells based on mitochondrial percentage reads, removal of doublets)
2. **Normalization**+integration+clustering + celltype marker
3. **Differential abundance** with MiloR
4. **Postclustering analysis**: Differential gene expression (for celltype marker), gene set enrichment, differential gene expression between 2 conditions, gene set enrichment.
5. **Cell-cell interaction**: MultinicheNet
6. **Trajectory inference**: Pseudotime with palantir, RNA velocity, CellRank
7. **Gene regulatory Networks**: SCENIC
