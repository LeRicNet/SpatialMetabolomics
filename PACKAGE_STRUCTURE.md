# SpatialMetabolics Package Structure

This document describes the complete file structure of the SpatialMetabolics R package.

## Root Directory Files

- `DESCRIPTION` - Package metadata, dependencies, and description
- `NAMESPACE` - Export and import directives (auto-generated from roxygen2)
- `LICENSE` - Package license (Artistic-2.0)
- `README.md` - Package overview and quick start guide
- `NEWS.md` - Version history and changelog
- `.Rbuildignore` - Files to exclude from package build
- `.gitignore` - Git version control exclusions
- `Makefile` - Development automation commands
- `SpatialMetabolics.Rproj` - RStudio project file
- `_pkgdown.yml` - Configuration for package website
- `CONTRIBUTING.md` - Contribution guidelines

## R/ Directory (Package Functions)

### Core Class and Methods
- `AllClasses.R` - S4 class definition for SpatialMetabolic
- `AllGenerics.R` - Generic method definitions
- `methods.R` - Accessor method implementations

### Data Import and Creation
- `createSpatialMetabolic.R` - Constructor functions

### Metabolic Pathways
- `metabolicPathways.R` - Curated gene sets for metabolic pathways

### Analysis Functions
- `calculateMetabolicScores.R` - Pathway activity scoring
- `spatialGradients.R` - Spatial pattern detection
- `compareStates.R` - Differential expression analysis

### Visualization
- `visualization.R` - Main plotting functions
- `visualization_extras.R` - Additional visualization utilities

### Utilities
- `utils.R` - Helper functions for data processing

## tests/ Directory (Unit Tests)

- `testthat.R` - Test runner
- `testthat/`
  - `test-SpatialMetabolic.R` - Core functionality tests
  - `test-pathways.R` - Metabolic pathway tests
  - `test-visualization.R` - Plot function tests

## inst/ Directory (Additional Files)

- `CITATION` - How to cite the package
- `script/`
  - `complete_analysis.R` - Full analysis workflow
  - `analyze_metabolic_spatial.R` - Basic usage example
- `extdata/`
  - `create_example_data.R` - Generate example datasets
- `shiny/`
  - `app.R` - Interactive Shiny application

## vignettes/ Directory (Long-form Documentation)

- `SpatialMetabolics.Rmd` - Main package vignette

## man/ Directory (Documentation)

Auto-generated from roxygen2 comments - contains .Rd files for all exported functions

## Key Features by File

### Data Structures
- **SpatialMetabolic class** (`AllClasses.R`): Extends SpatialExperiment with metabolic-specific slots

### Metabolic Pathways (`metabolicPathways.R`)
- OXPHOS complexes (I-V)
- Glycolysis pathway
- Mitochondrial biogenesis
- Mitochondrial dynamics
- ATP metabolism
- Small vessel disease signatures

### Analysis Methods
- **Scoring** (`calculateMetabolicScores.R`): Seurat, mean, ssGSEA methods
- **Spatial** (`spatialGradients.R`): Moran's I, Geary's C, variogram
- **Differential** (`compareStates.R`): Wilcoxon, t-test, permutation

### Visualization Types
- Spatial heatmaps
- Violin/box plots
- Volcano plots
- Pathway heatmaps
- Network plots
- QC plots
- Multi-panel figures

## Development Workflow

1. **Edit code** in R/ directory
2. **Document** with `make document`
3. **Test** with `make test`
4. **Check** with `make check`
5. **Build** with `make build`

## Installation

```r
# From GitHub
devtools::install_github("yourusername/SpatialMetabolics")

# From source
install.packages("SpatialMetabolics_1.0.0.tar.gz", repos = NULL)
```

## Main Dependencies

- **Bioconductor**: SpatialExperiment, SingleCellExperiment, BiocParallel
- **CRAN**: Seurat, ggplot2, viridis, patchwork
- **Suggests**: ComplexHeatmap, pheatmap, ggrepel

## Package Size

- **Source**: ~500 KB (excluding example data)
- **Installed**: ~2 MB
- **Example data**: ~10 MB (optional)

This structure follows Bioconductor package guidelines and provides a complete framework for spatial metabolic analysis.
