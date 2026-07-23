#==============================================================================
# Script: Figure5_plot_code.R
# Description: This script generates Figure 5 showing expression heatmaps of
#              photosynthesis-related genes across multiple tissues,
#              developmental stages, and environmental treatments for
#              representative Araliaceae species.
#==============================================================================


#==============================================================================
# 1. Generate Expression Heatmaps of Photosynthesis-related Genes
#==============================================================================


#------------------------------------------------------------------------------
# Set working directory
#------------------------------------------------------------------------------
setwd("D:/genomes/00figures_tables/Figure5/heatmap_photosynthesis240615/new_replace3")


#------------------------------------------------------------------------------
# Load required library
#------------------------------------------------------------------------------
library(pheatmap)


#------------------------------------------------------------------------------
# Example: Generate heatmap for A. continentalis
#------------------------------------------------------------------------------
# Input matrix:
#   - Rows: photosynthesis-related genes
#   - First column: module assignment
#   - Remaining columns: normalized expression values
#
# Genes are ordered according to module assignment.
# Samples are annotated by developmental stage, tissue, or treatment.

test <- read.table(
  "acofpkmraw.hotosynthesis.txt",
  stringsAsFactors = FALSE,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)


#------------------------------------------------------------------------------
# Define sample annotation
#------------------------------------------------------------------------------
fc <- factor(
  rep(
    c(
      "0h","12h","5d",
      "Shypocotyl","Sleaf","Spetiole",
      "Hhypocotyl","Hleaf","Hpetiole",
      "Lhypocotyl","Lleaf","Lpetiole",
      "Whypocotyl","Wleaf","Wpetiole"
    ),
    c(rep(1,15))
  ),
  levels = c(
    "0h","12h","5d",
    "Shypocotyl","Sleaf","Spetiole",
    "Hhypocotyl","Hleaf","Hpetiole",
    "Lhypocotyl","Lleaf","Lpetiole",
    "Whypocotyl","Wleaf","Wpetiole"
  )
)


#------------------------------------------------------------------------------
# Prepare data
#------------------------------------------------------------------------------
# Sort genes according to module assignment
test_sorted <- test[order(test$Module), ]

annotation_col <- data.frame(
  Environment = fc
)

rownames(annotation_col) <- colnames(test_sorted[-1])


#------------------------------------------------------------------------------
# Generate heatmap
#------------------------------------------------------------------------------
pdf(
  "acofpkmraw.hotosynthesis.pdf"
)

pheatmap(
  test_sorted[-1],
  color = colorRampPalette(
    c(
      "lightgray",
      "blue",
      "cadetblue3",
      "wheat",
      "hotpink",
      "red"
    )
  )(50),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 8,
  fontsize_row = 8,
  scale = "none",
  cellheight = 15,
  annotation_col = annotation_col,
  annotation_row = test_sorted[1],
  show_rownames = TRUE,
  border = FALSE
)

dev.off()



#==============================================================================
# 2. Generate Heatmaps for the Remaining Species
#==============================================================================
#
# The same procedure is repeated for the remaining Araliaceae species using
# their corresponding expression matrices and sample annotations:
#
#   - A. elata (Ael)
#   - Eleutherococcus senticosus (Ese)
#   - Daucus carota (Dca)
#   - Kalopanax septemlobus (Kse)
#   - Panax ginseng (Pgi)
#   - Panax quinquefolius (Pqu)
#   - Schefflera arboricola (Sar)
#   - Tetrapanax papyrifer (Tpa)
#
# Each species has its own sampling design (time points, tissues, and/or
# environmental treatments), while the visualization procedure remains
# identical:
#
#   1. Read the expression matrix.
#   2. Define sample annotations.
#   3. Sort genes by module assignment.
#   4. Generate a heatmap using pheatmap().
#   5. Export the figure as a PDF.
#
#==============================================================================


#==============================================================================
# End of script
#==============================================================================