
## LOAD PACKAGES ##

# Packages for data processing
library(haven)            # For reading in Stata format data (.dta)
library(tidyverse)        # For data manipulation and wrangling
library(labelled)         # For dealing with labels from Stata data
library(Hmisc)            # For weighted quantiles 
library(gt)               # For tables
library(gtsummary)        # For summary tables

# Packages for sequence analysis
library(TraMineR)         # For sequence analysis
library(seqHMM)           # For multi-channel sequence analysis
library(RColorBrewer)     # For colour palettes 
library(colorspace)       # For adjusting colour palettes
library(WeightedCluster)  # For cluster functions and statistics
library(grid)             # For gridplot
library(fpc)              # For bootstrapped cluster allocation

# Packages for sensitivity analysis
library(WeightIt)         # For inverse probability of attrition weighting
library(cobalt)           # For weighting diagnostics
library(Rmisc)            # For SF-12 confidence intervals

# Packages for multinomial logistic regression 
library(parsnip)          # For tidy model interfaces
library(broom)            # For tidy model output
library(broom.helpers)    # For additional tidying of outputs
library(yardstick)        # For model performance
library(brglm2)           # For bias reduction
library(ggh4x)            # To extend ggplots with additional features (facet_wrap colours)
library(cvms)             # For confusion matrices
library(patchwork)        # For combining plots together
library(flextable)        # For clean regression table outputs
