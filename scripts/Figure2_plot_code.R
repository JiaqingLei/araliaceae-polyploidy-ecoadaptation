#==============================================================================
# Script: Figure2_plot_code.R
# Description: This script generates Figure 2 visualizations:
#              1. Area charts showing gene number changes of paired ancestral
#                 chromosomes across species
#              2. Heatmaps showing GO category distributions across species
#==============================================================================


#==============================================================================
# 1. Area Chart Showing Gene Number Changes of Paired Ancestral Chromosomes
#==============================================================================


# Set working directory
setwd("D:\\genomes\\00figures_tables\\Figure2")


#------------------------------------------------------------------------------
# Load required libraries
#------------------------------------------------------------------------------
library(ggplot2)
library(data.table)
library(tidyr)
library(cowplot)   # Merge figures


#------------------------------------------------------------------------------
# Load gene number summary data
#------------------------------------------------------------------------------
# Input file contains the number of genes assigned to different ancestral
# chromosome categories across species.
#
# Type:
#   A01_Core          : conserved regions
#   A02_Dispensable   : dispensable regions
#   A03_Unique        : species-specific regions
#   B01_contraction   : contracted regions
#   B02_expansion     : expanded regions
#   C01_r_contraction : relative contraction regions
#   C02_r_expansion   : relative expansion regions
#
# Num2:
#   Number of genes assigned to each category

test <- read.table(
  "summary4_sort_negative_2.txt",
  header = TRUE,
  sep = "\t"
)



#------------------------------------------------------------------------------
# Select chromosome categories used for visualization
#------------------------------------------------------------------------------
# Retain ancestral chromosome categories related to gene conservation,
# expansion, and contraction patterns

dta1 <- subset(
  test,
  test[, "Type"] == "A01_Core" |
  test[, "Type"] == "A02_Dispensable" |
  test[, "Type"] == "A03_Unique" |
  test[, "Type"] == "B01_contraction" |
  test[, "Type"] == "B02_expansion" |
  test[, "Type"] == "C01_r_contraction" |
  test[, "Type"] == "C02_r_expansion"
)



#------------------------------------------------------------------------------
# Generate area chart
#------------------------------------------------------------------------------
# Each line represents a paired ancestral chromosome.
# The y-axis indicates gene number.
# Different colors represent paired ancestral chromosomes.
# Facets show different chromosome evolutionary categories.

plot.part1 <- ggplot(
  dta1,
  aes(
    x = Spp,
    y = Num2,
    fill = Chr,
    group = Chr,
    class = Mark
  )
) +
  geom_point(
    aes(color = Chr),
    shape = 15,
    size = 2
  ) +
  geom_segment(
    aes(
      x = Spp,
      xend = Spp,
      y = 0,
      yend = Num2,
      color = Chr
    )
  ) +
  geom_line(
    aes(color = Chr),
    linewidth = 1
  ) +
  geom_area(
    aes(fill = Chr),
    alpha = 0.3
  ) +
  facet_grid(
    Mark ~ Type
  ) +
  theme_classic() +
  scale_color_manual(
    values = c(
      "#1c9e77","#d95f02",
      "#1c9e77","#d95f02",
      "#1c9e77","#d95f02",
      "#1c9e77","#d95f02"
    )
  ) +
  scale_fill_manual(
    values = c(
      "#1c9e77","#d95f02",
      "#1c9e77","#d95f02",
      "#1c9e77","#d95f02",
      "#1c9e77","#d95f02"
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    col = "black"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 0.5
    )
  ) +
  ylab("Number of genes") +
  xlab("Species")


# Save figure
ggsave(
  file = "summary4_sort_negative_part1.pdf",
  plot.part1,
  width = 25,
  height = 20
)

ggsave(
  file = "summary4_sort_negative_part1.png",
  plot.part1,
  width = 25,
  height = 20
)


plot.part1




#==============================================================================
# 2. Heatmap Showing GO Category Distribution Across Species
#==============================================================================


# Set working directory
setwd("D:\\genomes\\00figures_tables\\Figure2\\new_figure2b")


#------------------------------------------------------------------------------
# Load required libraries
#------------------------------------------------------------------------------
library(pheatmap)
library(cowplot)
library(ggplotify)



#------------------------------------------------------------------------------
# Load GO enrichment matrix
#------------------------------------------------------------------------------
# Rows represent GO categories.
# Columns represent species.
# Values indicate the abundance of genes associated with each GO category.

data <- read.table(
  "A0.merge.go2_A01",
  stringsAsFactors = FALSE,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)



#------------------------------------------------------------------------------
# Remove GO categories without variation
#------------------------------------------------------------------------------
# Retain GO categories showing differences among species

test_A0_A01 <- data[
  apply(
    data,
    1,
    function(x) sd(x) != 0
  ),
]



#------------------------------------------------------------------------------
# Generate GO heatmap
#------------------------------------------------------------------------------
# Rows represent GO categories.
# Columns represent species.
# Row clustering groups similar functional patterns.

plot.A0_A01 <- as.ggplot(
  pheatmap(
    test_A0_A01,
    colorRampPalette(
      c(
        "#008792",
        "#00a6ac",
        "#78cdd1",
        "#dec674",
        "#c1a173",
        "#da765b"
      )
    )(6),
    cluster_rows = TRUE,
    treeheight_row = 0,
    cluster_cols = FALSE,
    scale = "row",
    fontsize = 6,
    show_rownames = FALSE,
    show_colnames = TRUE,
    angle_col = 90,
    legend = FALSE
  )
)


plot.A0_A01



#==============================================================================
# End of script
#==============================================================================

