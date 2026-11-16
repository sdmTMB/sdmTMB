# Keep the prepared data, parameters, map, and random list as the single model
# specification. Every caller that changes data builds a new objective.
# `adreport` optionally names the only ADREPORTs to register, so that
# `sdreport()` skips standard errors the caller won't use.
# `backend` is ignored; it is accepted because released sdmTMBextra passes it.
make_sdmTMB_adfun <- function(data, parameters, map, random = NULL,
                              backend = NULL, profile = NULL,
                              silent = TRUE, adreport = NULL, ...) {
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
  prepared <- rtmb_prepare(data)
  rtmb_validate(data, prepared, parameters, random, ...)
  objective <- rtmb_make_objective(prepared, adreport)
  RTMB::MakeADFun(objective, parameters = parameters, map = map,
    random = random, profile = profile, silent = silent, ...)
}
