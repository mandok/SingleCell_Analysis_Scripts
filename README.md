# SingleCell_Analysis_Scripts

This repository contains R/Python scripts for the analysis of scRNAseq and spatial transcriptomics


## scRNA-seq analysis
Inside the scRNAseq folder, the scripts covering the following analysis steps:
1. **Quality control** (removal of low quality cells based on mitochondrial percentage reads, removal of doublets)
2. **Normalization**+integration+clustering + celltype marker
3. **Differential abundance** with MiloR
4. **Postclustering analysis**: Differential gene expression (for celltype marker), gene set enrichment, differential gene expression between 2 conditions, gene set enrichment.
5. **Cell-cell interaction**: MultinicheNet
6. **Trajectory inference**: Pseudotime with palantir, RNA velocity, CellRank
7. **Gene regulatory Networks**: SCENIC


# Spatial Transcriptomics analysis
Inside the spatial_transcriptomics folder, the scripts cover the following analysis steps:
- Visium (spot resolution)
- Xenium (single-cell resolution)
- RESOLVE Biosciences (single-cell resolution)
