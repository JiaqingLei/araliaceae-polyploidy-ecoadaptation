# Analysis Protocol for Figure Reproduction

> **Associated Manuscript:** *Divergent genome fractionation and adaptation following polyploidy in the Araliaceae*

## Overview

This repository contains the R scripts used to reproduce the figures presented in the manuscript. To demonstrate the workflow while avoiding redistribution of complete research datasets, only example input datasets are provided. These example files illustrate the required input format and allow users to test each script. The complete datasets used to generate the publication figures are not included.

## Repository structure

```
araliaceae-polyploidy-ecoadaptation
├── example_data/
│    ├── Figure1/
│    ├── Figure2/
│    ├── Figure3/
│    ├── Figure4/
│    ├── Figure5/
├── scripts/
├── Protocol.md
└── README.md
```

## Software requirements

R 4.4.3

Required packages include ggplot2, tidyr, treemapify, cowplot, data.table, pheatmap, ggplotify, reshape2, dplyr, Rmisc, ComplexHeatmap, circlize and grid.

## Example data

Representative example datasets are provided in the example_data directory. These files demonstrate the expected input format, column organization and file naming convention required by each script. Users should replace the example files with their own processed data generated following the Materials and Methods section of the manuscript.

### 1. Figure1

**repreat_gene_other_final.txt**: Genome component composition matrix used for the percentage stacked bar chart, including repeat classes, gene components, and other genomic components for each species.
**length_pie.txt**: Input matrix used to generate species-specific pie charts. Each column corresponds to one species and contains values used for pie chart visualization.
**radius.txt**: Genome size values used to scale the radius of species-specific pie charts.
**5dup.txt**: Genomic duplication category matrix used for treemap visualization, including different duplication types across species.

### 2. Figure2

**summary4_sort_negative_2.txt**: Gene number matrix assigned to different ancestral chromosome categories across species.
**A0.merge.go2_A01**: GO category abundance matrix for ancestral chromosome-derived regions across species.

### 3. Figure3

**Ac.csv**: Photosynthetic rate measurements under different PPFD values used for fitting light-response curves.
**SLA.txt**: Specific leaf area measurements across species and experimental conditions.
**H_length.txt**: Phenotypic trait measurements related to hydraulic characteristics.
**P_length.txt**: Phenotypic trait measurements related to hydraulic characteristics.
**Chl_content.txt**: Chlorophyll content measurements across species and experimental conditions.
**00FvFm_14_herb-wood_4.txt**: Chlorophyll fluorescence measurements used for Fv/Fm analysis.

### 4. Figure4

**new_summary.final.merge_plot.txt**: Expression diversity and specificity matrix containing gene class and species information.
**Fig4C_Normal_light_leaf_FPKM.txt**: Gene expression matrix under normal light condition.
**Fig4C_Fluctuating_Shade_High_Low_light_leaf_DEG.txt**: DEG pattern matrix corresponding to the four variable light treatments.

### 5. Figure5

**acofpkmraw.hotosynthesis.txt**: Gene expression matrix of photosynthesis-related genes of Araliaceae species.
- Row names: photosynthesis-related genes
- First column: gene module assignment
- Remaining columns: normalized expression values across samples

## General workflow

- Generate processed summary tables from genome annotation, comparative genomics, transcriptome and other downstream analyses.
- Prepare input tables with the same structure as the example datasets.
- Place the processed files into the appropriate example_data/FigureX directory (or modify the input path in the script).
- Run the corresponding R script after updating the working directory if necessary.
- Preliminary draft of the corresponding figures are written to the output directory.
- Manually adjust and polish these draft figures to publication-quality standards.

## Reproducing figures

Each figure has an independent R script. Open the script, update the input path if required, replace the example dataset with the corresponding processed dataset, and execute the script using R or Rscript.

## Notes

The example datasets are intended solely to demonstrate the workflow and input format. They are not the complete datasets used in the manuscript. All plotting parameters are contained within the scripts.

## Reproducibility

The scripts were tested under R 4.4.3 on Windows and Linux. With properly formatted processed datasets, the scripts will reproduce preliminary draft of the corresponding figures.

## Citation

If you use these scripts or processed datasets, please cite our paper.
