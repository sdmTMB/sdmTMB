# Keep the prepared data, parameters, map, and random list as the single model
# specification. Every caller that changes data builds a new objective.
make_sdmTMB_adfun <- function(data, parameters, map, random = NULL,
                              backend = "tmb", profile = NULL,
                              silent = TRUE, ...) {
  backend <- match.arg(backend, c("tmb", "rtmb"))
  # The C++ template has no preferential-sampling likelihood; refuse rather
  # than fit the catch model alone.
  if (backend == "tmb" && has_preferential(data)) {
    cli::cli_abort("Preferential sampling requires backend = \"rtmb\"; set control = sdmTMBcontrol(backend = \"rtmb\").")
  }
  if (backend == "tmb") {
    obj <- TMB::MakeADFun(data = data, parameters = parameters, map = map,
      random = random, profile = profile, DLL = "sdmTMB", silent = silent, ...)
  } else {
    prepared <- rtmb_prepare(data)
    rtmb_validate(data, prepared, parameters, random, ...)
    objective <- rtmb_make_objective(prepared)
    obj <- RTMB::MakeADFun(objective, parameters = parameters, map = map,
      random = random, profile = profile, silent = silent, ...)
  }
  attr(obj, "sdmTMB_backend") <- backend
  obj
}

sdreport_sdmTMB <- function(obj, ...) {
  backend <- attr(obj, "sdmTMB_backend")
  if (is.null(backend) || backend == "tmb") return(TMB::sdreport(obj, ...))
  RTMB::sdreport(obj, ...)
}

# Fits saved before backend metadata existed used the C++ TMB model. Once the
# C++ model is removed, those fits will need a refit instead.
backend_sdmTMB <- function(object) {
  if (is.null(object$backend)) "tmb" else object$backend
}
