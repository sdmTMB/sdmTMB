# tinyVAST SAR model example
# --------------------------
# This script fits a simultaneous autoregressive (SAR) spatial model
# using tinyVAST's built-in salmon_returns example data.

library(tinyVAST)
library(igraph)

options("tinyVAST.verbose" = FALSE)

# Built-in example data
data(salmon_returns)

# Remove zeros for a lognormal response example
salmon_returns$Biomass_nozeros <- ifelse(
  salmon_returns$Biomass == 0,
  NA,
  salmon_returns$Biomass
)

Data <- na.omit(salmon_returns)

# SAR adjacency graph among regions
adjacency_graph <- make_graph(
  ~ Korea - Japan - M.I - WKam - EKam -
    WAK - SPen - Kod - CI - PWS -
    SEAK - NBC - SBC - WA
)

# AR(2) dynamics for each salmon species
dsem <- "
  sockeye -> sockeye, -1, lag1_sockeye
  sockeye -> sockeye, -2, lag2_sockeye

  pink -> pink, -1, lag1_pink
  pink -> pink, -2, lag2_pink

  chum -> chum, -1, lag1_chum
  chum -> chum, -2, lag2_chum
"

# Fit model with SAR spatial structure
fit_sar <- tinyVAST(
  formula = Biomass_nozeros ~ 0 + Species + Region,
  data = Data,
  spacetime_term = dsem,
  variable_column = "Species",
  time_column = "Year",
  space_columns = "Region",
  distribution_column = "Species",
  family = list(
    chum = lognormal(),
    pink = lognormal(),
    sockeye = lognormal()
  ),
  spatial_domain = adjacency_graph,
  control = tinyVASTcontrol(profile = "alpha_j")
)

# Summarize spatiotemporal parameters
summary(fit_sar, what = "spacetime_term")

# Model AIC
AIC(fit_sar)
