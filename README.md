# Analysis scripts for human blood sciMET data

[![DOI](https://zenodo.org/badge/1393414675.svg)](https://doi.org/10.5281/zenodo.23337538)

> **Disclaimer:** These scripts were written for manuscript-specific analyses and are provided for transparency and reproducibility. They are not intended to function as a standalone or production-ready software package.

## Contents

- `QC`: notebook for initial sciMET QC, using the PBMC dataset as an example.
- `analysis`: code for sciMET data preprocessing, including various data format conversion, doublet detection, and metagene profile computation etc.
- `matrix_clustering`: notebook for MethSCAn matrix computation, downstream dimensionality reduction and clustering analysis on VMR matrix.
- `genotype_profile`: genetic variant calling, processing and dimensionality reduction.
- `differential_methylation`: scripts for MethSCAn DMR identification and processing.
- `methylvi_integration`: scripts for running MethylVI data integration with the public single-cell DNA methylation dataset.
- `scIB_benchmark`: scripts for data preparation and computation of scIB benchmarker metrics.
- `replication_timing`: scripts for replication timing analysis.
