#==============================================================================
# Script: Figure4C_plot_code.R
# Description: This script generates Figure 4C showing:
#              1. FPKM expression heatmap for Normal light condition
#              2. DEG heatmaps for Fluctuating, Shade, High, and Low light conditions
#              All heatmaps are combined vertically for comparative visualization
#==============================================================================
rm(list = ls())

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(grid)

#==============================================================================
# 1. Define Global Constants
#==============================================================================
# Define column order once for all heatmaps (consistent across all conditions)
COLUMN_ORDER <- c(
  "Kse|Ara_sub1", "Kse|Ara_sub2", "Kse|Arb_sub1", "Kse|Arb_sub2",     # K. septemlobus
  "Ese|Ara_sub1", "Ese|Ara_sub2", "Ese|Arb_sub1", "Ese|Arb_sub2",     # E. senticosus
  "Tpa|Ara_sub1", "Tpa|Ara_sub2", "Tpa|Arb_sub1", "Tpa|Arb_sub2",     # T. papyrifer
  "Sar|Ara_sub1", "Sar|Ara_sub2", "Sar|Arb_sub1", "Sar|Arb_sub2",     # S. arboricola
  "Ael|Ara", "Ael|Arb",                                               # A. elata                                          # A. elata
  "Aco|Ara_sub1", "Aco|Ara_sub2", "Aco|Arb_sub1", "Aco|Arb_sub2",     # A. continentalis
  "Pqu|Ara_sub1", "Pqu|Ara_sub2", "Pqu|Arb_sub1", "Pqu|Arb_sub2",     # P. quinquefolius
  "Pgi|Ara_sub1", "Pgi|Ara_sub2", "Pgi|Arb_sub1", "Pgi|Arb_sub2"      # P. ginseng
)

# Define gene order for consistent visualization
GENE_ORDER <- c(
  "CRY1", "CRY2", 
  "PHOT1", "PHOT2",
  "PHYA", "PHYB", "PHYC", "PHYE",
  "PIF1", "PIF3", "PIF7", "PIF8",
  "HY5"
)

# Define column text colors based on experimental groups
# Red/Pink for Ara-derived lines, Blue/Lightblue for Arb-derived lines
COL_COLORS <- c(
  rep(c("red", "pink", "blue", "lightblue"), 4),
  c("red", "blue"),
  rep("red", 2),
  rep("blue", 2),
  rep(c("red", "pink", "blue", "lightblue"), 2)
)

# Color mapping for DEG: 0=no change, 1=down-regulated, 2=up-regulated
DEG_COLORS <- c("0" = "white", "1" = "blue", "2" = "red", "NA" = "lightgray")

#==============================================================================
# 2. Define Helper Functions
#==============================================================================
#' Process DEG data for a specific condition
#' @param deg_data Raw DEG data frame
#' @param condition_name Name of the condition column (e.g., "Fluctuation", "Shade")
#' @param gene_order Vector of gene names in desired order
#' @param column_order Vector of column names in desired order
#' @return A matrix with genes as rows and conditions as columns
process_deg_data <- function(deg_data, condition_name, gene_order, column_order) {
  # Aggregate values by Symbol and Idchr
  deg_group <- deg_data %>%
    select(Symbol, Idchr, all_of(condition_name)) %>%
    group_by(Symbol, Idchr) %>%
    summarise(!!sym(condition_name) := sum(!!sym(condition_name)), .groups = "drop") %>%
    arrange(Idchr)
  
  # Reshape to wide format
  deg_wide <- deg_group %>%
    pivot_wider(
      names_from = Idchr, 
      values_from = all_of(condition_name),
      values_fill = NA
    ) %>%
    select(Symbol, all_of(column_order))
  
  # Convert to data frame and reorder rows
  deg_df <- as.data.frame(deg_wide)
  deg_df <- deg_df[match(gene_order, deg_df$Symbol), ]
  
  # Remove rows that are all NA (genes not present in this condition)
  deg_df <- deg_df[!apply(is.na(deg_df[-1]), 1, all), ]
  
  # Set row names and remove Symbol column
  rownames(deg_df) <- deg_df$Symbol
  deg_df <- deg_df[-1]
  
  return(deg_df)
}

#' Create a DEG heatmap with standardized parameters
#' @param data_matrix Matrix of DEG values
#' @param output_file PDF file name for output
#' @param col_colors Column name color vector
#' @param deg_colors Color mapping for DEG values
#' @param column_split Factor for column splitting
#' @return A Heatmap object
create_deg_heatmap <- function(data_matrix, output_file, col_colors, deg_colors, column_split) {
  pdf(output_file, width = 10, height = 6)
  
  ht <- Heatmap(
    data_matrix,
    col = deg_colors,
    na_col = "lightgray",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 9, col = col_colors),
    row_names_side = "left",
    name = "DEG",
    column_split = column_split,
    column_title = NULL,
    width = unit(15, "cm"),
    height = unit(2, "cm"),
    column_gap = unit(1, "mm"),
    row_gap = unit(10, "mm"),
    border = TRUE,
    border_gp = gpar(col = "darkred", lwd = 1.5, lty = 1),
    heatmap_legend_param = list(
      at = c(1, 2),
      labels = c("down", "up")
    )
  )
  
  draw(ht)
  dev.off()
  
  return(ht)
}

#==============================================================================
# 3. Data Loading and Preprocessing - FPKM Expression Data
#==============================================================================
## Load FPKM expression data - Normal light condition
wlfpkm <- read.table("Fig4C_Normal_light_leaf_FPKM.txt", header = TRUE, sep = "\t")

## Aggregate FPKM values by Type and Idchr (summing replicates/isoforms)
wlfpkm_group <- wlfpkm %>%
  group_by(Type, Idchr) %>%
  summarise(FPKM_sum = sum(FPKM), .groups = "drop") %>%
  arrange(Idchr)

## Reshape data from long to wide format
wkfpkm_groupnew <- wlfpkm_group %>%
  select(Type, Idchr, FPKM_sum) %>%
  pivot_wider(
    names_from = Idchr,
    values_from = FPKM_sum,
    values_fill = NA
  ) %>%
  select(Type, all_of(COLUMN_ORDER))

## Convert to matrix format with Type as row names
wlfpkm_gndf <- as.data.frame(wkfpkm_groupnew)
rownames(wlfpkm_gndf) <- wlfpkm_gndf$Type
wlfpkm_gndf <- wlfpkm_gndf[-1]

#==============================================================================
# 4. Generate FPKM Expression Heatmap
#==============================================================================
## Create column split factor (each column as separate group)
col_group_factor <- factor(1:ncol(wlfpkm_gndf))

## Convert to matrix and apply log2 transformation
wlfpkm_mat <- as.matrix(log2(wlfpkm_gndf + 1))

## Build FPKM heatmap with white-to-red color gradient
pdf("Fig4C_Normal_light_leaf_FPKM_heatmap.pdf", width = 10, height = 6)
htW <- Heatmap(
  wlfpkm_mat,
  col = colorRamp2(c(0, 10), c("white", "red")),
  na_col = "darkgray",
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9, col = COL_COLORS),
  row_names_side = "left",
  column_names_side = "top",
  name = "log2(FPKM+1)",
  column_split = col_group_factor,
  column_title = NULL,
  width = unit(15, "cm"),
  height = unit(2, "cm"),
  column_gap = unit(1, "mm"),
  row_gap = unit(10, "mm"),
  border = TRUE,
  border_gp = gpar(col = "darkred", lwd = 1.5, lty = 1)
)
draw(htW)
dev.off()

#==============================================================================
# 5. Load DEG Data and Process All Conditions
#==============================================================================
## Load DEG data containing the other four light conditions
deg <- read.table("Fig4C_Fluctuating_Shade_High_Low_light_leaf_DEG.txt", 
                  header = TRUE, sep = "\t")

## Define all conditions to process
conditions <- c("Fluctuation", "Shade", "High", "Low")
condition_names <- c("Fluctuating", "Shade", "High", "Low")
output_files <- c(
  "Fig4C_Fluctuating_light_DEG_heatmap.pdf",
  "Fig4C_Shade_light_DEG_heatmap.pdf",
  "Fig4C_High_light_DEG_heatmap.pdf",
  "Fig4C_Low_light_DEG_heatmap.pdf"
)

## Process each condition and store heatmap objects
ht_list <- list(htW)  # Start with FPKM heatmap

for (i in seq_along(conditions)) {
  # Process DEG data for this condition
  deg_df <- process_deg_data(
    deg_data = deg,
    condition_name = conditions[i],
    gene_order = GENE_ORDER,
    column_order = COLUMN_ORDER
  )
  
  # Convert to matrix
  deg_mat <- as.matrix(deg_df)
  
  # Create heatmap
  ht <- create_deg_heatmap(
    data_matrix = deg_mat,
    output_file = output_files[i],
    col_colors = COL_COLORS,
    deg_colors = DEG_COLORS,
    column_split = col_group_factor
  )
  
  # Store for combined figure
  ht_list <- c(ht_list, list(ht))
}

#==============================================================================
# 6. Combine All Heatmaps into a Single Figure
#==============================================================================
## Vertically stack all five heatmaps (FPKM + 4 DEG conditions)
pdf("Fig4C_Combined_Normal_light_FPKM_Fluctuating_Shade_High_Low_light_DEG_heatmap.pdf", 
    width = 10, height = 8)
combined_heatmap <- Reduce("%v%", ht_list)
draw(combined_heatmap, row_gap = unit(10, "mm"))
dev.off()

cat("\n=== Script completed successfully ===\n")
#==============================================================================
# End of script
#==============================================================================
