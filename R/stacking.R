#' Perform stacking with log scores on `sdmTMB_cv()` output
#'
#' @description
#' This approach is described in Yao et al. (2018) \doi{10.1214/17-BA1091}. The
#' general method minimizes (or maximizes) some quantity across models. For
#' simple models with normal error, this may be the root mean squared error
#' (RMSE), but other approaches include the log score. We adopt the latter here,
#' where log scores are used to generate the stacking of predictive
#' distributions
#'
#' @param model_list A list of models fit with [sdmTMB_cv()] to generate
#'   estimates of predictive densities. You will want to set the seed
#'   to the same value before fitting each model or manually construct
#'   the fold IDs so that they are the same across models.
#' @param include_folds An optional numeric vector specifying which folds to
#'   include in the calculations. For example, if 5 folds are used for k-fold
#'   cross validation, and the first 4 are needed to generate these weights,
#'   `include_folds = 1:4`.
#' @importFrom stats optim runif
#' @export
#' @return
#' A vector of model weights.
#'
#' @references
#' Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. 2018. Using Stacking to
#' Average Bayesian Predictive Distributions (with Discussion). Bayesian Analysis
#' 13(3): 917–1007. International Society for Bayesian Analysis.
#' \doi{10.1214/17-BA1091}
#'
#'
#' @examples
#' \donttest{
#' # Set parallel processing if desired. See 'Details' in ?sdmTMB_cv
#'
#' # Depth as quadratic:
#' set.seed(1)
#' m_cv_1 <- sdmTMB_cv(
#'   density ~ 0 + depth_scaled + depth_scaled2,
#'   data = pcod_2011, mesh = pcod_mesh_2011,
#'   family = tweedie(link = "log"), k_folds = 2
#' )
#' # Depth as linear:
#' set.seed(1)
#' m_cv_2 <- sdmTMB_cv(
#'   density ~ 0 + depth_scaled,
#'   data = pcod_2011, mesh = pcod_mesh_2011,
#'   family = tweedie(link = "log"), k_folds = 2
#' )
#'
#' # Only an intercept:
#' set.seed(1)
#' m_cv_3 <- sdmTMB_cv(
#'   density ~ 1,
#'   data = pcod_2011, mesh = pcod_mesh_2011,
#'   family = tweedie(link = "log"), k_folds = 2
#' )
#'
#' models <- list(m_cv_1, m_cv_2, m_cv_3)
#' weights <- sdmTMB_stacking(models)
#' weights
#' }
sdmTMB_stacking <- function(model_list, include_folds = NULL) {
  n_models <- length(model_list)
  if (is.null(include_folds)) {
    n_folds <- max(model_list[[1]]$data$cv_fold)
    include_folds <- seq_len(n_folds)
  }

  flds <- lapply(model_list, function(x) x$data$cv_fold)
  flds_chk <- lapply(flds[seq(2, length(flds))], function(x) identical(flds[[1]], x))
  flds_same <- all(unlist(flds_chk))
  if (!flds_same) {
    cli_abort("Not all models used the same folds.")
  }

  # the only quantity we need is the log likelihood
  X <- matrix(0,
    nrow = nrow(model_list[[1]]$data),
    ncol = n_models
  )
  for (i in seq_len(n_models)) X[, i] <- model_list[[i]]$data$cv_loglik

  # filter out only records in folds of interest
  X <- X[model_list[[1]]$data$cv_fold %in% include_folds, , drop = FALSE]

  # exponentiate to convert to density in normal space, eqn 5 in Yao et al. (2018)
  X <- exp(X)

  # find model weights that maximize the likelihood
  tot_ll <- function(p, X) {
    # softmax / normalize to compositions
    z <- matrix(exp(p) / sum(exp(p)), ncol = 1L)
    # see equation below eq 5 in Yao et al. sum over data points of log product of
    # weights and probabilities
    -sum(log(X %*% z)) # neg X because optim will minimize NLL
  }
  # use optim to find weights that maximize total
  o <- optim(par = rep(0, n_models), fn = tot_ll, X = X)
  # softmax / normalize to compositions
  weights <- exp(o$par) / sum(exp(o$par))
  weights
}
