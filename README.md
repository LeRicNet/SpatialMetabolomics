# SpatialMetabolics

## Spatial Metabolic Analysis for Visium HD Brain Transcriptomics

SpatialMetabolics is an R/Bioconductor package designed for comprehensive analysis of ATP regulation, mitochondrial function, and metabolic states in spatial transcriptomics data. It is optimized for Visium HD brain tissue analysis with seamless Seurat integration.

## Features

- **Metabolic Pathway Analysis**: Pre-defined gene sets for oxidative phosphorylation, glycolysis, mitochondrial biogenesis/dynamics, and ATP metabolism
- **Spatial Pattern Detection**: Identify metabolic gradients and spatially variable pathways using Moran's I and other spatial statistics
- **Differential Analysis**: Compare metabolic states between conditions with spatial-aware statistical methods
- **Publication-Ready Visualization**: Create high-quality figures with minimal code
- **Seurat Integration**: Seamless conversion between Seurat and SpatialMetabolic objects
- **Memory Efficient**: Optimized for large Visium HD datasets

## Installation

### Prerequisites

```r
# Install Bioconductor if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c(
    "SpatialExperiment",
    "SingleCellExperiment", 
    "SummarizedExperiment",
    "BiocParallel",
    "sparseMatrixStats",
    "ComplexHeatmap"
))

# Install CRAN dependencies
install.packages(c(
    "Seurat",
    "ggplot2", 
    "viridis",
    "patchwork",
    "pheatmap",
    "ggrepel",
    "cowplot"
))
```

### Install SpatialMetabolics

```r
# Install from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("yourusername/SpatialMetabolics")
```

## Quick Start

```r
library(SpatialMetabolics)
library(Seurat)

# Load Visium HD data
seurat_wt <- Load10X_Spatial("data/WT_sample/", bin.size = 8)
seurat_svd <- Load10X_Spatial("data/SVD_sample/", bin.size = 8)

# Add metadata
seurat_wt$condition <- "WT"
seurat_wt$sample <- "WT1"
seurat_svd$condition <- "SVD"
seurat_svd$sample <- "SVD1"

# Create SpatialMetabolic object
spm <- createSpatialMetabolic(
    list(seurat_wt, seurat_svd),
    condition_col = "condition",
    sample_col = "sample"
)

# Normalize and calculate metabolic scores
spm <- normalizeSpatial(spm)
spm <- calculateMetabolicScores(spm, 
    pathways = c("oxphos", "glycolysis", "atp"))

# Detect spatial patterns
spm <- detectMetabolicGradients(spm)

# Compare conditions
results <- compareMetabolicStates(spm,
    ref_group = "WT",
    test_group = "SVD")

# Visualize
p <- plotSpatialMetabolicScore(spm, pathway = "OXPHOS")
print(p)
```

## Main Functions

### Data Import and Processing
- `createSpatialMetabolic()`: Create SpatialMetabolic object from Seurat objects
- `normalizeSpatial()`: Normalize spatial transcriptomics data
- `findVariableFeatures()`: Identify highly variable genes
- `runPCA()`: Principal component analysis

### Metabolic Analysis
- `calculateMetabolicScores()`: Calculate pathway activity scores
- `detectMetabolicGradients()`: Find spatially variable metabolic patterns
- `compareMetabolicStates()`: Differential metabolic analysis
- `findMetabolicNeighborhoods()`: Identify metabolically distinct regions

### Pathway Gene Sets
- `getOXPHOSGenes()`: Oxidative phosphorylation complex genes
- `getGlycolysisGenes()`: Glycolysis pathway genes
- `getMitoBiogenesisGenes()`: Mitochondrial biogenesis regulators
- `getMitoDynamicsGenes()`: Mitochondrial fission/fusion genes
- `getATPGenes()`: ATP synthesis and consumption genes
- `getSVDGenes()`: Small vessel disease signature genes

### Visualization
- `plotSpatialMetabolicScore()`: Spatial heatmaps of pathway activity
- `plotMetabolicComparison()`: Compare pathways between conditions
- `plotSpatialGradients()`: Top spatially variable features
- `plotMetabolicHeatmap()`: Sample-level pathway summaries
- `plotVolcano()`: Volcano plots for differential expression

## Example Analysis: Small Vessel Disease

```r
# Full analysis pipeline
library(SpatialMetabolics)
library(patchwork)

# Load and process data
spm <- createSpatialMetabolic(seurat_list)
spm <- preprocessSpatial(spm, n_variable_features = 3000)

# Calculate all metabolic pathways
pathways <- getAllMetabolicPathways("mouse")
spm <- calculateMetabolicScores(spm, pathways = pathways)

# Spatial analysis
spm <- detectMetabolicGradients(spm)

# Statistical comparison
de_pathways <- compareMetabolicStates(spm, 
    ref_group = "WT", test_group = "SVD")

# Create publication figure
p1 <- plotSpatialMetabolicScore(spm, "OXPHOS", 
    samples = c("WT1", "SVD1"))
p2 <- plotMetabolicComparison(spm, "OXPHOS")
p3 <- plotSpatialGradients(spm, top_n = 10)
p4 <- plotVolcano(de_pathways)

figure <- (p1 | p2) / (p3 | p4) +
    plot_annotation(tag_levels = 'A')

ggsave("metabolic_analysis.pdf", figure, 
    width = 10, height = 10)
```

## Advanced Features

### Custom Pathway Analysis
```r
# Define custom pathways
my_pathways <- list(
    MyPathway1 = c("Gene1", "Gene2", "Gene3"),
    MyPathway2 = c("Gene4", "Gene5", "Gene6")
)

spm <- calculateMetabolicScores(spm, pathways = my_pathways)
```

### Metabolic Neighborhoods
```r
# Find metabolically distinct regions
spm <- findMetabolicNeighborhoods(spm, k = 5)

# Visualize neighborhoods
plotSpatialFeatures(spm, 
    features = "metabolic_neighborhood",
    color_scale = "Set3")
```

### Integration with Other Tools
```r
# Convert to Seurat for additional analysis
seurat_obj <- toSeurat(spm, add_metabolic_scores = TRUE)

# Use Seurat functions
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj)
```

## Troubleshooting

### Memory Issues
For large Visium HD datasets:
```r
# Use parallel processing
library(BiocParallel)
spm <- calculateMetabolicScores(spm,
    BPPARAM = MulticoreParam(workers = 4))

# Process in chunks
samples <- unique(colData(spm)$sample_id)
results <- lapply(samples, function(s) {
    spm_sub <- subsetSpatial(spm, samples = s)
    compareMetabolicStates(spm_sub)
})
```

### Missing Genes
Check pathway coverage:
```r
# List available pathways and gene counts
pathway_info <- listMetabolicPathways("mouse")
print(pathway_info)

# Check which genes are present
pathways <- getAllMetabolicPathways("mouse")
genes_present <- lapply(pathways, function(x) {
    intersect(x, rownames(spm))
})
```

## Citation

If you use SpatialMetabolics in your research, please cite:

```
Your Name et al. (2024). SpatialMetabolics: Spatial metabolic analysis 
for Visium HD brain transcriptomics. R package version 1.0.0.
https://github.com/yourusername/SpatialMetabolics
```

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## License

Artistic-2.0

## Contact

- Bug reports: https://github.com/yourusername/SpatialMetabolics/issues
- Email: your.email@example.com

## Acknowledgments

This package builds upon the excellent work of:
- Bioconductor SpatialExperiment developers
- Seurat development team
- 10x Genomics Visium HD platform
