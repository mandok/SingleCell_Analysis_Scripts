# SingleCell_Analysis_Scripts

This repository contains R/Python scripts for the analysis of scRNAseq and spatial transcriptomics


## scRNA-seq analysis
Inside the scRNAseq folder, the scripts covering the following analysis steps using [Seurat V5](https://satijalab.org/seurat/):
1. **Quality control** (removal of low quality cells based on mitochondrial percentage reads, removal of doublets with [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
2. **Normalization**: 
3. **Integration**
4. **Clustering**
5. **Differential abundance** with [MiloR](https://bioconductor.org/packages/release/bioc/vignettes/miloR/inst/doc/milo_demo.html)
6. **Postclustering analysis**: Differential gene expression (for celltype marker), gene set enrichment, differential gene expression between 2 conditions, gene set enrichment.
7. **Cell-cell interaction**: MultinicheNet
8. **Trajectory inference**: Pseudotime with palantir, RNA velocity, CellRank
9. **Gene regulatory Networks**: SCENIC


# Spatial Transcriptomics analysis
Inside the spatial_transcriptomics folder, the scripts cover the following analysis steps:
- Visium (spot resolution)
- Xenium (single-cell resolution)
- RESOLVE Biosciences (single-cell resolution)
