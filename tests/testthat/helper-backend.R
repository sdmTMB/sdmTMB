backend <- Sys.getenv("SDMTMB_TEST_BACKEND")
if (backend %in% c("tmb", "rtmb")) {
  options(sdmTMB.backend = backend)
}

# The censored binomial, beta-binomial, and NB families are RTMB only, so their
# tests use RTMB even when `SDMTMB_TEST_BACKEND` selects TMB.
local_rtmb_backend <- function(env = parent.frame()) {
  withr::local_options(sdmTMB.backend = "rtmb", .local_envir = env)
}
