#==============================================================================
# Script: Figure3_plot_code.R
# Description: This script generates Figure 3 phenotypic and physiological trait
#              visualizations:
#              1. Light-response curves showing photosynthetic responses under
#                 different light intensities
#              2. Boxplots showing variation in phenotypic traits and
#                 chlorophyll content among species and conditions
#              3. Boxplots showing chlorophyll fluorescence parameter Fv/Fm
#==============================================================================


#==============================================================================
# 1. Light-response Curves
#==============================================================================


#------------------------------------------------------------------------------
# Load required libraries
#------------------------------------------------------------------------------
library(reshape2)
library(ggplot2)



#------------------------------------------------------------------------------
# Load photosynthetic response data
#------------------------------------------------------------------------------
# Input file:
#   Ac.csv
#
# The dataset contains photosynthetic rate measurements under different
# photosynthetic photon flux density (PPFD) levels for different samples.
#
# Columns:
#   PPFD   : Photosynthetic photon flux density
#   Sample : Sample identity
#   Value  : Photosynthetic rate

lrc <- read.csv("Ac.csv")



#------------------------------------------------------------------------------
# Reshape data from wide format to long format
#------------------------------------------------------------------------------
# Convert samples into a single column format for ggplot visualization.

lrc <- melt(
  lrc,
  id = "PPFD"
)

colnames(lrc) <- c(
  "PPFD",
  "Sample",
  "value"
)



#------------------------------------------------------------------------------
# Fit non-linear light-response model
#------------------------------------------------------------------------------
# Photosynthetic response curves are fitted using a non-linear model:
#
# Pn = alpha*((1 - beta*x)/(1 + gamma*x))*x - Rd
#
# Parameters:
#   alpha : apparent quantum efficiency
#   beta  : photoinhibition coefficient
#   gamma : curvature parameter
#   Rd    : dark respiration rate

p <- ggplot(
  data = lrc,
  aes(
    x = PPFD,
    y = value,
    group = Sample,
    color = Sample,
    shape = Sample
  )
)


p1 <- p +
  geom_smooth(
    method = "nls",
    formula = y ~ alpha*((1 - beta*x)/(1 + gamma*x))*x - Rd,
    se = FALSE,
    method.args = list(
      start = c(
        alpha = 0.07,
        beta = 0.00005,
        gamma = 0.004,
        Rd = 0.2
      )
    ),
    linewidth = 3
  ) +
  labs(
    y = expression(paste("Photosynthetic rate")),
    x = expression(paste("PAR"))
  )



#------------------------------------------------------------------------------
# Format light-response curves
#------------------------------------------------------------------------------

p1 +
  scale_x_continuous(
    breaks = seq(0, 1000, by = 200)
  ) +
  scale_y_continuous(
    limits = c(-2, 8),
    breaks = seq(-2, 8, 2)
  ) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "#CAB2D6",
      "#B2DF8A",
      "#FAA43A",
      "#1F78B4"
    )
  ) +
  scale_fill_manual(
    values = c(
      "#CAB2D6",
      "#B2DF8A",
      "#FAA43A",
      "#1F78B4"
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  )



#==============================================================================
# 2. Boxplots of Phenotypic Traits and Chlorophyll Content
#==============================================================================


#------------------------------------------------------------------------------
# Load required libraries
#------------------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Rmisc)



#------------------------------------------------------------------------------
# Define color palette for experimental conditions
#------------------------------------------------------------------------------

mypal <- c(
  "#6666FF",
  "#FF6633",
  "aquamarine2",
  "darkgoldenrod2",
  "chocolate1",
  "#FFC619"
)



#------------------------------------------------------------------------------
# Function:
# Generate boxplots showing phenotypic trait variation
#------------------------------------------------------------------------------
# Each boxplot compares phenotypic trait values among species and conditions.
# Individual points represent biological replicates.


#------------------------------------------------------------------------------
# 2.1 Specific leaf area (SLA)
#------------------------------------------------------------------------------
# Input file:
#   SLA.txt
#
# SLA represents a morphological trait related to leaf structure.

dta <- read.table(
  "SLA.txt",
  header = TRUE,
  sep = "\t",
  row.names = NULL,
  check.names = FALSE
)

dta1 <- gather(
  dta,
  key = "Code",
  value = "num",
  3:ncol(dta)
)



plot1 <- ggplot(
  dta1,
  aes(
    x = Species,
    y = num,
    fill = Condition,
    color = Condition
  )
) +
  geom_boxplot(
    width = 0.5,
    linewidth = 0.2,
    color = "black",
    outlier.alpha = 0,
    position = position_dodge(width = 0.6)
  ) +
  geom_jitter(
    alpha = 0.5,
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.6
    )
  ) +
  scale_color_manual(
    values = mypal
  ) +
  scale_fill_manual(
    values = mypal
  ) +
  ylab("SLA") +
  xlab("Species") +
  theme(
    legend.position = c(0.95,0.7),
    axis.text.x = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5
    ),
    panel.background = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = "transparent"
    )
  )

plot1



#------------------------------------------------------------------------------
# 2.2 H_length trait
#------------------------------------------------------------------------------
# Input file:
#   H_length.txt

dta <- read.table(
  "H_length.txt",
  header = TRUE,
  sep = "\t",
  row.names = NULL,
  check.names = FALSE
)

dta1 <- gather(
  dta,
  key = "Code",
  value = "num",
  3:ncol(dta)
)



plot2 <- ggplot(
  dta1,
  aes(
    x = Species,
    y = num,
    fill = Condition,
    color = Condition
  )
) +
  geom_boxplot(
    width = 0.5,
    linewidth = 0.2,
    color = "black",
    outlier.alpha = 0,
    position = position_dodge(width = 0.6)
  ) +
  geom_jitter(
    alpha = 0.5,
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.6
    )
  ) +
  scale_color_manual(values = mypal) +
  scale_fill_manual(values = mypal) +
  ylab("H_length") +
  xlab("Species")

plot2



#------------------------------------------------------------------------------
# 2.3 P_length trait
#------------------------------------------------------------------------------
# Input file:
#   P_length.txt


dta <- read.table(
  "P_length.txt",
  header = TRUE,
  sep = "\t",
  row.names = NULL,
  check.names = FALSE
)

dta1 <- gather(
  dta,
  key = "Code",
  value = "num",
  3:ncol(dta)
)



plot3 <- ggplot(
  dta1,
  aes(
    x = Species,
    y = num,
    fill = Condition,
    color = Condition
  )
) +
  geom_boxplot(
    width = 0.5,
    linewidth = 0.2,
    color = "black",
    outlier.alpha = 0,
    position = position_dodge(width = 0.6)
  ) +
  geom_jitter(
    alpha = 0.5,
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.6
    )
  ) +
  scale_color_manual(values = mypal) +
  scale_fill_manual(values = mypal) +
  ylab("P_length") +
  xlab("Species")

plot3



#------------------------------------------------------------------------------
# 2.4 Chlorophyll content
#------------------------------------------------------------------------------
# Input file:
#   Chl_content.txt


dta <- read.table(
  "Chl_content.txt",
  header = TRUE,
  sep = "\t",
  row.names = NULL,
  check.names = FALSE
)

dta1 <- gather(
  dta,
  key = "Code",
  value = "num",
  3:ncol(dta)
)



plot4 <- ggplot(
  dta1,
  aes(
    x = Species,
    y = num,
    fill = Condition,
    color = Condition
  )
) +
  geom_boxplot(
    width = 0.5,
    linewidth = 0.2,
    color = "black",
    outlier.alpha = 0,
    position = position_dodge(width = 0.6)
  ) +
  geom_jitter(
    alpha = 0.5,
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.6
    )
  ) +
  scale_color_manual(values = mypal) +
  scale_fill_manual(values = mypal) +
  ylab("Chl_content") +
  xlab("Species")

plot4



#------------------------------------------------------------------------------
# Combine phenotypic trait plots
#------------------------------------------------------------------------------

plot <- plot_grid(
  plot1,
  plot2,
  plot3,
  plot4,
  ncol = 1,
  nrow = 4
)


ggsave(
  file = "Figure.pdf",
  plot,
  width = 10,
  height = 15
)



#==============================================================================
# 3. Boxplot of Chlorophyll Fluorescence Parameter Fv/Fm
#==============================================================================


#------------------------------------------------------------------------------
# Load Fv/Fm dataset
#------------------------------------------------------------------------------
# Input file:
#   00FvFm_14_herb-wood_4.txt
#
# Fv/Fm represents the maximum quantum efficiency of photosystem II (PSII).

setwd("G:/FvFm")


dta <- read.table(
  "00FvFm_14_herb-wood_4.txt",
  header = TRUE,
  sep = "\t",
  row.names = NULL,
  check.names = FALSE
)


#------------------------------------------------------------------------------
# Reshape Fv/Fm data
#------------------------------------------------------------------------------

dta1 <- gather(
  dta,
  key = "Code",
  value = "num",
  3:ncol(dta)
)



#------------------------------------------------------------------------------
# Generate Fv/Fm boxplot
#------------------------------------------------------------------------------

plot1 <- ggplot(
  dta1,
  aes(
    x = Species,
    y = num,
    fill = Condition,
    color = Condition
  )
) +
  geom_boxplot(
    width = 0.5,
    linewidth = 0.2,
    color = "black",
    outlier.alpha = 0,
    position = position_dodge(width = 0.6)
  ) +
  geom_jitter(
    alpha = 0.5,
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.6
    )
  ) +
  scale_color_manual(values = mypal) +
  scale_fill_manual(values = mypal) +
  ylab("Fv/Fm value") +
  xlab("Species")


plot1



#------------------------------------------------------------------------------
# Save Fv/Fm figure
#------------------------------------------------------------------------------

ggsave(
  file = "00FvFm_14_herb-wood_4.pdf",
  plot1,
  width = 12,
  height = 10
)



#==============================================================================
# End of script
#==============================================================================