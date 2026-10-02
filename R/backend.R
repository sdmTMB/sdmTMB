# Keep the prepared data, parameters, map, and random list as the single model
# specification. Every caller that changes data builds a new objective.
# `adreport` optionally names the only ADREPORTs to register (RTMB backend
# only), so that `sdreport()` skips standard errors the caller won't use.
make_sdmTMB_adfun <- function(data, parameters, map, random = NULL,
                              backend = "tmb", profile = NULL,
                              silent = TRUE, adreport = NULL, ...) {
  backend <- match.arg(backend, c("tmb", "rtmb"))
  # The C++ template has no preferential-sampling likelihood; refuse rather
  # than fit the catch model alone.
  if (backend == "tmb" && has_preferential(data)) {
    cli::cli_abort("Preferential sampling requires backend = \"rtmb\"; set control = sdmTMBcontrol(backend = \"rtmb\").")
  }
  # Parameters of removed experimental epsilon models; fits saved before the
  # removal still carry them (mapped off unless the option was used).
  removed <- intersect(names(parameters),
    c("b_epsilon", "epsilon_re", "ln_epsilon_re_sigma"))
  if (length(removed)) {
    if (!all(vapply(map[removed], function(x) length(x) && all(is.na(x)), logical(1)))) {
      cli_abort(c("This model uses a removed experimental `epsilon_model` option.",
        "i" = "Refit it with the installed sdmTMB version."))
    }
    parameters[removed] <- NULL
    map[removed] <- NULL
    if (!is.null(random)) random <- setdiff(random, removed)
  }
  # Weighted-average vectors are one per column (x and y for COG); fits
  # saved before this carry a plain vector
  if (!is.null(data$proj_vector)) data$proj_vector <- as.matrix(data$proj_vector)
  if (backend == "tmb") {
    obj <- TMB::MakeADFun(data = data, parameters = parameters, map = map,
      random = random, profile = profile, DLL = "sdmTMB", silent = silent, ...)
  } else {
    prepared <- rtmb_prepare(data)
    rtmb_validate(data, prepared, parameters, random, ...)
    objective <- rtmb_make_objective(prepared, adreport)
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
