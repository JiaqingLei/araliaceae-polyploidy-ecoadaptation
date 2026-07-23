#==============================================================================
# Script: Figure1_plot_code.R
# Description: This script generates Figure 1 genome composition visualizations:
#              1. Percentage stacked bar chart showing repeat and gene composition
#              2. Species-specific pie charts showing genome size
#              3. Treemap plots showing different genomic duplication categories
#==============================================================================

#==============================================================================
# 1. Percentage Stacked Bar Chart Showing Repeat and Gene Composition
#==============================================================================

# Set working directory
setwd("D:\\genomes\\00figures_tables\\Figure1")


#------------------------------------------------------------------------------
# Load required libraries
#------------------------------------------------------------------------------
library(ggplot2)
library(tidyr)


#------------------------------------------------------------------------------
# Load repeat and gene composition data
#------------------------------------------------------------------------------
# Input file contains genome component proportions
# (repeat classes, genes, and other genomic components) for each species
dta <- read.table(
  "repreat_gene_other_final.txt",
  header = TRUE,
  sep = "\t"
)


#------------------------------------------------------------------------------
# Reshape data from wide format to long format
#------------------------------------------------------------------------------
# Convert species-specific columns into a single column for ggplot visualization
dta1 <- gather(
  dta,
  key = "spp",
  value = "num",
  2:21
)


#------------------------------------------------------------------------------
# Define colors for genome components
#------------------------------------------------------------------------------
mypal2 <- c(
  "cadetblue2",
  "cornflowerblue",
  "dodgerblue1",
  "blue",
  "grey",
  "darkorange4",
  "chocolate1",
  "coral",
  "darkgoldenrod2"
)


#------------------------------------------------------------------------------
# Generate percentage stacked bar chart
#------------------------------------------------------------------------------
# Each bar represents one species;
# different colors indicate different genome components

p3 <- ggplot(
  dta1,
  aes(
    x = spp,
    y = num,
    fill = Class
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.8,
    colour = "black"
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2)
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  scale_fill_manual(
    values = mypal2
  )


# Display plot
p3



#==============================================================================
# 2. Species-specific Pie Charts Showing Genome Size
#==============================================================================

# Set working directory
setwd("D:\\genomes\\00figures_tables\\Figure1")


#------------------------------------------------------------------------------
# Load genome size data
#------------------------------------------------------------------------------
# dat_radius: genome size values used to scale pie chart radius
# dat_length: values used to generate individual pie charts

dat_length <- read.table(
  "length_pie.txt",
  header = TRUE,
  sep = "\t"
)

dat_radius <- read.table(
  "radius.txt",
  header = TRUE,
  sep = "\t"
)


#------------------------------------------------------------------------------
# Generate species-specific pie charts
#------------------------------------------------------------------------------
# Each pie chart represents one species.
# Pie radius is proportional to genome size.
# Labels and borders are removed for cleaner visualization.

pdf(
  file = "length_pie.pdf",
  width = 6,
  height = 10
)

par(
  mfcol = c(20, 1),
  mar = c(0.1, 0.1, 0.1, 0.1),
  oma = c(0.1, 0.1, 0.1, 0.1)
)


for (f in 2:21) {
  
  pie(
    dat_length[, f],
    radius = dat_radius[, f],
    labels = NA,
    border = NA,
    col = "darkorange1"
  )
}


dev.off()


#==============================================================================
# 3. Treemap Visualization of Genomic Duplication Categories
#==============================================================================


#------------------------------------------------------------------------------
# Load required packages
#------------------------------------------------------------------------------
# Install treemapify if necessary:
# devtools::install_github("wilkox/treemapify")

library(ggplot2)
library(treemapify)
library(cowplot)


# Set working directory
setwd("D:\\genomes\\00figures_tables\\Figure1")


#------------------------------------------------------------------------------
# Load treemap input data
#------------------------------------------------------------------------------
# Input file contains different genomic duplication categories
# for each species

dat_5dup <- read.table(
  "5dup.txt",
  header = TRUE,
  sep = "\t"
)



#------------------------------------------------------------------------------
# Define function to generate treemap plots
#------------------------------------------------------------------------------
# Each treemap represents genomic duplication composition for one species.
# Area size corresponds to the abundance of each duplication category.
# Colors indicate different duplication types.

create_treemap <- function(column_index) {
  
  ggplot(
    data = dat_5dup,
    aes(
      area = dat_5dup[, column_index],
      fill = Parameters,
      label = Parameters
    )
  ) +
    geom_treemap() +
    geom_treemap_text() +
    scale_fill_manual(
      values = c(
        "dodgerblue1",
        "grey",
        "chartreuse3",
        "purple",
        "pink",
        "darkorange1"
      )
    ) +
    scale_x_continuous(
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    )
}


#------------------------------------------------------------------------------
# Generate treemap plots for each species
#------------------------------------------------------------------------------
# Columns 2-21 correspond to different species

treemap_list <- list()

for (i in 2:21) {
  
  treemap_list[[i]] <- create_treemap(i)
  
}


#------------------------------------------------------------------------------
# Combine treemap plots
#------------------------------------------------------------------------------
# Merge plots vertically for figure assembly

plot_grid(
  plotlist = treemap_list[2:8],
  nrow = 7,
  ncol = 1
)

plot_grid(
  plotlist = treemap_list[9:15],
  nrow = 7,
  ncol = 1
)

plot_grid(
  plotlist = treemap_list[16:21],
  nrow = 6,
  ncol = 1
)



#==============================================================================
# End of script
#==============================================================================
