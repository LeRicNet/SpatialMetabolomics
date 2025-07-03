# SpatialMetabolics 1.0.0

## New Features

* Initial release of SpatialMetabolics package
* S4 class `SpatialMetabolic` extending `SpatialExperiment` for metabolic analysis
* Comprehensive metabolic pathway gene sets for mouse and human
  - Oxidative phosphorylation (OXPHOS) complexes I-V
  - Glycolysis pathway
  - Mitochondrial biogenesis and dynamics
  - ATP synthesis and consumption
  - Small vessel disease signatures

## Analysis Methods

* `calculateMetabolicScores()`: Pathway activity scoring with Seurat method
* `detectMetabolicGradients()`: Spatial pattern detection using Moran's I
* `compareMetabolicStates()`: Differential metabolic analysis between conditions
* `findMetabolicNeighborhoods()`: Identify metabolically distinct spatial regions

## Visualization

* `plotSpatialMetabolicScore()`: Spatial heatmaps with customizable color scales
* `plotMetabolicComparison()`: Violin/box plots for condition comparisons
* `plotSpatialGradients()`: Top spatially variable features
* `plotMetabolicHeatmap()`: Sample-level pathway activity heatmaps
* `plotVolcano()`: Volcano plots with automatic labeling
* Additional plots: network, correlation, QC, and summary reports

## Data Processing

* `createSpatialMetabolic()`: Convert from Seurat objects
* `normalizeSpatial()`: Multiple normalization methods
* `findVariableFeatures()`: Feature selection
* `runPCA()` and `runUMAP()`: Dimension reduction

## Utilities

* Seamless integration with Seurat workflows
* Memory-efficient processing for Visium HD data
* Parallel processing support via BiocParallel
* Comprehensive documentation and vignette
