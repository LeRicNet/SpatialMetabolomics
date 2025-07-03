# SpatialMetabolics Installation Guide

## System Requirements

- **R version**: ≥ 4.1.0
- **Bioconductor version**: ≥ 3.14
- **Operating System**: Windows, macOS, or Linux
- **Memory**: ≥ 8 GB RAM recommended for Visium HD data
- **Disk Space**: ~100 MB for package and dependencies

## Quick Installation

### From Bioconductor (when available)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SpatialMetabolics")
```

### From GitHub (development version)

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("yourusername/SpatialMetabolics")
```

## Full Installation with Dependencies

### 1. Install Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = "3.18")  # Use latest version
```

### 2. Install Core Dependencies

```r
# Bioconductor packages
BiocManager::install(c(
    "SpatialExperiment",
    "SingleCellExperiment",
    "SummarizedExperiment",
    "S4Vectors",
    "BiocGenerics",
    "BiocParallel",
    "sparseMatrixStats",
    "DelayedArray"
))

# CRAN packages
install.packages(c(
    "Seurat",
    "SeuratObject",
    "ggplot2",
    "viridis",
    "patchwork",
    "scales",
    "cowplot",
    "Matrix"
))
```

### 3. Install Optional Dependencies

```r
# For enhanced visualization
install.packages(c(
    "pheatmap",
    "ggrepel",
    "RColorBrewer",
    "ggridges",
    "plotly"
))

BiocManager::install("ComplexHeatmap")

# For the Shiny app
install.packages(c(
    "shiny",
    "shinydashboard",
    "DT"
))

# For network plots
install.packages(c(
    "igraph",
    "ggraph"
))

# For additional functionality
install.packages(c(
    "uwot",      # UMAP
    "edgeR",     # TMM normalization
    "testthat",  # Running tests
    "knitr",     # Building vignettes
    "rmarkdown", # Building vignettes
    "BiocStyle"  # Vignette styling
))

BiocManager::install("scran")  # scran normalization
```

### 4. Install SpatialMetabolics

```r
# From source file
install.packages("SpatialMetabolics_1.0.0.tar.gz", repos = NULL, type = "source")

# Or from GitHub
devtools::install_github("yourusername/SpatialMetabolics", 
                        build_vignettes = TRUE)
```

## Troubleshooting Installation

### Common Issues

#### 1. Dependency conflicts

```r
# Check for conflicts
BiocManager::valid()

# Update all packages
BiocManager::install(ask = FALSE)
```

#### 2. Compilation errors (macOS/Linux)

Ensure you have development tools installed:
- **macOS**: Install Xcode Command Line Tools
  ```bash
  xcode-select --install
  ```
- **Linux**: Install R development packages
  ```bash
  # Ubuntu/Debian
  sudo apt-get install r-base-dev
  
  # Fedora/RHEL
  sudo yum install R-devel
  ```

#### 3. Memory issues during installation

```r
# Increase memory limit (Windows)
memory.limit(size = 8000)

# Install packages one at a time
deps <- c("SpatialExperiment", "Seurat", "ggplot2")
for (pkg in deps) {
    BiocManager::install(pkg)
    gc()  # Garbage collection
}
```

#### 4. Seurat installation issues

```r
# Install Seurat dependencies first
install.packages(c("spatstat.explore", "spatstat.geom"))
install.packages("Seurat")
```

## Verifying Installation

```r
# Load the package
library(SpatialMetabolics)

# Check version
packageVersion("SpatialMetabolics")

# Run example
example(createSpatialMetabolic)

# Run tests
testthat::test_package("SpatialMetabolics")
```

## Docker Installation

For a reproducible environment, use our Docker image:

```dockerfile
FROM bioconductor/bioconductor_docker:RELEASE_3_18

RUN R -e "BiocManager::install('yourusername/SpatialMetabolics')"
```

Build and run:
```bash
docker build -t spatialmetabolics .
docker run -it spatialmetabolics R
```

## Platform-Specific Notes

### Windows
- Rtools is required for building from source
- Download from: https://cran.r-project.org/bin/windows/Rtools/

### macOS
- May need to install gfortran for some dependencies
- Silicon (M1/M2) Macs: Use native R build for best performance

### Linux
- Additional system libraries may be required:
  ```bash
  # Ubuntu/Debian
  sudo apt-get install \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev
  ```

## Post-Installation Setup

### 1. Set up BiocParallel for your system

```r
library(BiocParallel)

# Check available cores
BiocParallel::bpparam()

# Set up for your system
# For Unix/macOS
register(MulticoreParam(workers = 4))

# For Windows
register(SnowParam(workers = 4))
```

### 2. Configure memory settings for large datasets

```r
# Increase memory limit (Windows)
memory.limit(size = 16000)

# Set DelayedArray block size
DelayedArray::setAutoBlockSize(1e9)  # 1GB blocks
```

### 3. Test with example data

```r
# Create small example dataset
library(SpatialMetabolics)

# Generate test data
spe <- SpatialExperiment:::.create_random_spe()
spm <- as(spe, "SpatialMetabolic")

# Run basic analysis
spm <- normalizeSpatial(spm)
pathways <- list(test = rownames(spm)[1:10])
spm <- calculateMetabolicScores(spm, pathways)

# Create a plot
plotSpatialMetabolicScore(spm, "test")
```

## Getting Help

- **Documentation**: `?SpatialMetabolics`
- **Vignette**: `vignette("SpatialMetabolics")`
- **GitHub Issues**: https://github.com/yourusername/SpatialMetabolics/issues
- **Bioconductor Support**: https://support.bioconductor.org/

## Citing SpatialMetabolics

```r
citation("SpatialMetabolics")
```
