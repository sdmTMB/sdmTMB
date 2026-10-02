backend <- Sys.getenv("SDMTMB_TEST_BACKEND")
if (backend %in% c("tmb", "rtmb")) {
  options(sdmTMB.backend = backend)
}
