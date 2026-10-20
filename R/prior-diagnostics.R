#' Inspect custom priors
#'
#' `get_prior_parameters()` lists the parameters that custom priors can refer
#' to: the elements of `par` and `theta` passed to the `custom` and
#' `custom_log_jacobian` functions in [sdmTMBpriors()], with their values. Use
#' it with a fitted model or one set up with `do_fit = FALSE`.
#'
#' `get_prior_densities()` evaluates the custom prior functions at the
#' estimated parameters, including random effects at their conditional modes,
#' and returns the log density of each term. Like the objective, it evaluates
#' `custom_log_jacobian` only if `bayesian = TRUE`.
#'
#' @details
#' `par` holds the raw parameters, with the names and shapes used internally.
#' Matrix columns are model components (e.g., the two parts of a delta model);
#' the second component's main effects are in `b_j2`. `theta` holds
#' natural-scale transformations of `par`, such as `phi = exp(ln_phi)` and
#' `range = sqrt(8) / kappa`.
#'
#' Elements of `par` that are fixed, by `map` or because the model doesn't
#' use them, are left out, since a prior on them would have no effect. All elements of `theta`
#' are listed, including transformations of parameters this model doesn't use
#' (e.g., `tweedie_p` in a model without a Tweedie family). Check that the
#' `par` elements underlying a `theta` element are listed before putting a
#' prior on it. Standard deviations of random fields the model doesn't
#' include (e.g., `sigma_O` without a spatial field) are zero.
#'
#' These names and shapes are internal and could change in future versions of
#' sdmTMB. Check them with `get_prior_parameters()` when writing custom priors.
#'
#' @param object An [sdmTMB()] model fit with the RTMB backend. For
#'   `get_prior_densities()`, it must have custom priors.
#' @param random Include the random effects (e.g., random field values)? There
#'   are often many of them.
#'
#' @return
#' `get_prior_parameters()`: a data frame with one row per element:
#' * `expression`: how to refer to the element in a custom prior function;
#' * `list`, `name`: `"par"` or `"theta"`, and the element's name;
#' * `label`: the coefficient name, where available;
#' * `value`: the estimate (or starting value if the model isn't fit);
#' * `status`: for `par`, `"estimated"`, `"random"`, or `"shared"` (estimated
#'   jointly with other elements via `map`); `"derived"` for `theta`.
#'
#' `get_prior_densities()`: a data frame with one row per term added to the
#' objective: `type` (`"density"`, or `"log_jacobian"` if `bayesian = TRUE`),
#' `term` (the returned name, or its position), and `log_density`.
#' @seealso [sdmTMBpriors()]
#' @export
#' @examples
#' fit <- sdmTMB(density ~ depth_scaled, data = pcod_2011,
#'   family = tweedie(), mesh = pcod_mesh_2011,
#'   control = sdmTMBcontrol(backend = "rtmb"),
#'   priors = sdmTMBpriors(custom = function(par, theta) {
#'     c(depth = RTMB::dnorm(par$b_j[2], 0, 1, log = TRUE))
#'   }))
#' get_prior_parameters(fit)
#' get_prior_densities(fit)
get_prior_parameters <- function(object, random = FALSE) {
  assert_that(inherits(object, "sdmTMB"))
  assert_that(is.logical(random), length(random) == 1L)
  if (backend_sdmTMB(object) != "rtmb") {
    cli_abort("Custom priors need a model fit with `sdmTMBcontrol(backend = \"rtmb\")`.")
  }
  prepared <- rtmb_prepare(object$tmb_data)
  par <- prior_parlist(object)
  par_rows <- flatten_prior_elements(par, "par")
  status <- rep("estimated", nrow(par_rows))
  status[par_rows$name %in% object$tmb_random] <- "random"
  for (n in intersect(names(par), names(object$tmb_map))) {
    map <- as.integer(object$tmb_map[[n]])
    i <- par_rows$name == n
    status[i][!is.na(map) & map %in% map[duplicated(map)]] <- "shared"
    status[i][is.na(map)] <- "fixed"
  }
  par_rows$status <- status
  theta_rows <- flatten_prior_elements(rtmb_transform(par, prepared), "theta")
  theta_rows$status <- rep("derived", nrow(theta_rows))
  out <- rbind(par_rows, theta_rows)
  out$label <- NA_character_
  for (m in seq_len(prepared$n_m)) {
    b <- out$list == "par" & out$name == c("b_j", "b_j2")[m]
    out$label[b] <- colnames(object$tmb_data$X_ij[[m]])
  }
  drop <- c("fixed", if (!random) "random")
  out <- out[!out$status %in% drop, c("expression", "list", "name", "label",
    "value", "status")]
  rownames(out) <- NULL
  out
}

#' @rdname get_prior_parameters
#' @export
get_prior_densities <- function(object) {
  assert_that(inherits(object, "sdmTMB"))
  custom <- object$tmb_data$priors_custom
  if (is.null(custom)) cli_abort("This model has no custom priors.")
  prepared <- rtmb_prepare(object$tmb_data)
  par <- prior_parlist(object)
  theta <- rtmb_transform(par, prepared)
  used <- c("density", if (prepared$priors$stan) "log_jacobian")
  rows <- lapply(used, function(which) {
    x <- rtmb_custom_terms(par, theta, custom, which)
    if (is.null(x)) return(NULL)
    term <- names(x)
    if (is.null(term)) term <- rep("", length(x))
    term[term == ""] <- as.character(seq_along(x))[term == ""]
    data.frame(type = which, term = term, log_density = unname(as.numeric(x)))
  })
  do.call(rbind, rows)
}

# Parameters as saved with the fit, or the starting values with
# `do_fit = FALSE`.
prior_parlist <- function(object) {
  if (!is.null(object$parlist)) object$parlist else object$tmb_obj$env$parList()
}

# One row per element of a (nested) list of arrays, with an `expression` for
# indexing it, e.g. `par$ln_kappa[1, 2]` or `theta$H[[1]][2, 1]`.
flatten_prior_elements <- function(x, prefix) {
  rows <- lapply(names(x), function(n) {
    flatten_prior_element(x[[n]], paste0(prefix, "$", n), n)
  })
  out <- do.call(rbind, rows)
  out$list <- rep(prefix, nrow(out))
  out
}

flatten_prior_element <- function(x, expression, name) {
  if (!length(x)) return(NULL)
  if (is.list(x)) {
    nested <- if (is.null(names(x))) {
      paste0(expression, "[[", seq_along(x), "]]")
    } else {
      paste0(expression, "$", names(x))
    }
    rows <- lapply(seq_along(x), function(i) {
      flatten_prior_element(x[[i]], nested[[i]], name)
    })
    return(do.call(rbind, rows))
  }
  x <- as.array(unclass(x))
  d <- dim(x)
  index <- if (length(d) == 1L) {
    as.character(seq_along(x))
  } else {
    apply(arrayInd(seq_along(x), d), 1L, paste, collapse = ", ")
  }
  data.frame(
    expression = paste0(rep(expression, length(x)), "[", index, "]"),
    name = rep(name, length(x)),
    value = as.numeric(x)
  )
}
