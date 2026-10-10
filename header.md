# sdmTMB

> Spatial and spatiotemporal GLMMs with TMB

[![CRAN
version](https://www.r-pkg.org/badges/version/sdmTMB)](https://cran.r-project.org/package=sdmTMB)
[![Documentation](https://img.shields.io/badge/documentation-sdmTMB-orange.svg?colorB=E91E63)](https://sdmTMB.github.io/sdmTMB/)
[![R-CMD-check](https://github.com/sdmTMB/sdmTMB/workflows/R-CMD-check/badge.svg)](https://github.com/sdmTMB/sdmTMB/actions)
[![Codecov test
coverage](https://codecov.io/gh/sdmTMB/sdmTMB/branch/main/graph/badge.svg)](https://app.codecov.io/gh/sdmTMB/sdmTMB?branch=main)
[![downloads](https://cranlogs.r-pkg.org/badges/sdmTMB)](https://cranlogs.r-pkg.org/)

sdmTMB is an R package that fits spatial and spatiotemporal GLMMs
(generalized linear mixed effects models) using Template Model Builder
([TMB](https://github.com/kaskr/adcomp) via
[RTMB](https://github.com/kaskr/RTMB)),
[fmesher](https://github.com/inlabru-org/fmesher), and Gaussian Markov
random fields. It’s designed to feel familiar to users of glm(), lme4,
mgcv, or glmmTMB. One common application is for species distribution
models (SDMs), but it works for any data with spatial coordinates. See
the [documentation site](https://sdmTMB.github.io/sdmTMB/) and the
[published paper](https://doi.org/10.18637/jss.v115.i02).

## Table of contents
