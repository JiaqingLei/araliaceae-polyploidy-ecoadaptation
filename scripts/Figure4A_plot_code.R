#==============================================================================
# Script: Figure4A_plot_code.R
# Description: This script generates Figure 4A showing the relationship between
#              expression diversity (Hj) and expression specificity (Sj)
#              across gene classes for representative species.
#==============================================================================


#==============================================================================
# 1. Diversity–Specificity Scatter Plots
#==============================================================================


#------------------------------------------------------------------------------
# Set working directory
#------------------------------------------------------------------------------
setwd("G:/NC/tpm")


#------------------------------------------------------------------------------
# Load required libraries
#------------------------------------------------------------------------------
library(ggplot2)
library(cowplot)



#------------------------------------------------------------------------------
# Load input data
#------------------------------------------------------------------------------
# Input file:
#   new_summary.final.merge_plot.txt
#
# The dataset contains expression diversity (Hj), expression specificity (Sj),
# gene class information, and species identity.

file <- list.files(pattern = "new_summary.final.merge_plot.txt")
n <- length(file)
file

dat<- read.table(file[1],header=T,sep="\t")

#------------------------------------------------------------------------------
# Generate diversity–specificity plot for A. elata
#------------------------------------------------------------------------------
# Each point represents one gene.
# x-axis: expression diversity (Hj)
# y-axis: expression specificity (Sj)
# Colors and symbols indicate different gene classes.
# The grey line represents the linear regression fit with 95% confidence interval.

data1 <- subset(dat, dat[,5] == "Ael")

plot_Ael <- ggplot(
  data1,
  aes(
    x = data1[,1],
    y = data1[,2],
    shape = class,
    color = class
  )
) +
  geom_point(
    size = 2,
    stroke = 1
  ) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1.5,
    fill = "grey",
    alpha = 0.2
  ) +
  scale_shape_manual(
    values = c(1,0,0,0,2,2,2)
  ) +
  scale_color_manual(
    values = c(
      "purple",
      "cadetblue3",
      "#7185EB",
      "blue",
      "red",
      "#FF8C00",
      "darkgoldenrod1"
    )
  ) +
  scale_x_continuous(
    limits = c(10,14),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    limits = c(0,0.65),
    expand = c(0,0)
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = "transparent"
    )
  ) +
  ggtitle("Ael") +
  xlab("Diversity (Hj)") +
  ylab("Specificity (Sj)")

plot_Ael


# Save figure
ggsave(
  file = paste0("0_", file[1], "_Ael.pdf"),
  plot_Ael,
  width = 5,
  height = 5
)



#------------------------------------------------------------------------------
# Generate diversity–specificity plot for P. ginseng
#------------------------------------------------------------------------------
# Visualization settings are identical to those used for A. elata.

data1 <- subset(dat, dat[,5] == "Pgi")

plot_Pgi <- ggplot(
  data1,
  aes(
    x = data1[,1],
    y = data1[,2],
    shape = class,
    color = class
  )
) +
  geom_point(
    size = 2,
    stroke = 1
  ) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1.5,
    fill = "grey",
    alpha = 0.2
  ) +
  scale_shape_manual(
    values = c(1,0,0,2,2,2)
  ) +
  scale_color_manual(
    values = c(
      "purple",
      "#7185EB",
      "blue",
      "red",
      "#FF8C00",
      "darkgoldenrod1"
    )
  ) +
  scale_x_continuous(
    limits = c(10,14),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    limits = c(0,0.65),
    expand = c(0,0)
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = "transparent"
    )
  ) +
  ggtitle("Pgi") +
  xlab("Diversity (Hj)") +
  ylab("Specificity (Sj)")

plot_Pgi


# Save figure
ggsave(
  file = paste0("0_", file[1], "_Pgi.pdf"),
  plot_Pgi,
  width = 5,
  height = 5
)



#==============================================================================
# End of script
#==============================================================================



