# Keep the prepared data, parameters, map, and random list as the single model
# specification. Every caller that changes data builds a new objective.
make_sdmTMB_adfun <- function(data, parameters, map, random = NULL,
                              backend = "tmb", profile = NULL,
                              silent = TRUE, ...) {
  backend <- match.arg(backend, c("tmb", "rtmb"))
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
