qres_tweedie <- function(object, y, mu, ...) {
  theta <- get_pars(object)
  p <- stats::plogis(theta[["thetaf"]]) + 1
  dispersion <- exp(theta[["ln_phi"]])

  u <- fishMod::pTweedie(q = y, p = p, mu = mu, phi = dispersion)
  if (p > 1 && p < 2) {
    u[y == 0] <- stats::runif(sum(y == 0), min = 0, max = u[y == 0])
  }
  stats::qnorm(u)
}

qres_binomial <- function(object, y, mu, .n = NULL) {
  # p <- object$family$linkinv(mu) # robust binomial in link space!
  p <- mu
  if (is.null(.n)) .n <- rep(1, length(y))
  mu <- .n * mu
  a <- stats::pbinom(y - 1, .n, p)
  b <- stats::pbinom(y, .n, p)
  u <- stats::runif(n = length(y), min = pmin(a, b), max = pmax(a, b))
  stats::qnorm(u)
}

# Beta-binomial CDF: P(X <= q) for X ~ BetaBinomial(size, alpha, beta)
# Replaces extraDistr::pbbinom to avoid orphaned package dependency
pbbinom <- function(q, size, alpha, beta) {
  # Vectorize over all inputs
  len <- max(length(q), length(size), length(alpha), length(beta))
  q <- rep_len(q, len)
  size <- rep_len(size, len)
  alpha <- rep_len(alpha, len)
  beta <- rep_len(beta, len)

  result <- numeric(len)

  for (i in seq_len(len)) {
    if (is.na(q[i]) || is.na(size[i]) || is.na(alpha[i]) || is.na(beta[i])) {
      result[i] <- NA
      next
    }

    # Handle edge cases
    if (q[i] < 0) {
      result[i] <- 0
      next
    }
    if (q[i] >= size[i]) {
      result[i] <- 1
      next
    }

    # Sum PMF from 0 to floor(q)
    # PMF: choose(n, k) * beta(k + alpha, n - k + beta) / beta(alpha, beta)
    # Using log scale for numerical stability
    q_int <- floor(q[i])
    n <- size[i]
    a <- alpha[i]
    b <- beta[i]
    log_beta_ab <- lbeta(a, b)
    cdf <- 0

    for (k in 0:q_int) {
      log_pmf <- lchoose(n, k) + lbeta(k + a, n - k + b) - log_beta_ab
      cdf <- cdf + exp(log_pmf)
    }

    result[i] <- cdf
  }

  result
}

qres_betabinomial <- function(object, y, mu, .n = NULL) {
  # Extract dispersion parameter
  theta <- get_pars(object)
  phi <- exp(theta[["ln_phi"]])
  if (is_delta(object)) phi <- phi[2]

  # TMB parameterization: alpha = mu * phi, beta = (1 - mu) * phi
  # where mu is already on probability scale from linkinv
  alpha <- mu * phi
  beta <- (1 - mu) * phi

  if (is.null(.n)) .n <- rep(1, length(y))

  a <- pbbinom(pmax(0, y - 1), size = .n, alpha = alpha, beta = beta)
  b <- pbbinom(y, size = .n, alpha = alpha, beta = beta)
  u <- stats::runif(n = length(y), min = pmin(a, b), max = pmax(a, b))
  stats::qnorm(u)
}

qres_nbinom2 <- function(object, y, mu, ...) {
  theta <- get_pars(object)
  phi <- exp(theta[["ln_phi"]])
  if (is_delta(object)) phi <- phi[2]
  a <- stats::pnbinom(y - 1, size = phi, mu = mu)
  b <- stats::pnbinom(y, size = phi, mu = mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}

# from glmmTMB:
rnbinom1 <- function(n, mu, phi) {
  # var = mu*(1+phi) = mu*(1+mu/k) -> k = mu/phi
  stats::rnbinom(n, mu = mu, size = mu / phi)
}
dnbinom1 <- function(x, mu, phi) {
  stats::dnbinom(x, mu = mu, size = mu / phi)
}
pnbinom1 <- function(q, mu, phi) {
  stats::pnbinom(q, mu = mu, size = mu / phi)
}
qnbinom1 <- function(p, mu, phi) {
  stats::qnbinom(p, mu = mu, size = mu / phi)
}

qres_nbinom1 <- function(object, y, mu, ...) {
  theta <- get_pars(object)
  phi <- exp(theta[["ln_phi"]])
  if (is_delta(object)) phi <- phi[2]
  a <- pnbinom1(y - 1, phi = phi, mu = mu)
  b <- pnbinom1(y, phi = phi, mu = mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}

ptruncated_nbinom2 <- function(q, mu, phi){
  p <- (stats::pnbinom(q, size = phi, mu = mu) - stats::pnbinom(0, size = phi, mu = mu))/(1 - stats::pnbinom(0, size = phi, mu = mu))
}
ptruncated_nbinom1 <- function(q, mu, phi){
  p <- (pnbinom1(q, mu, phi) - pnbinom1(0, mu, phi))/(1 - pnbinom1(0, mu, phi))
}

qres_truncated_nbinom2 <- function(object, y, mu, ...) {
  theta <- get_pars(object)
  phi <- exp(theta[["ln_phi"]])
  if (is_delta(object)) phi <- phi[2]
  a <- ptruncated_nbinom2(y - 1, mu = mu, phi = phi)
  b <- ptruncated_nbinom2(y, mu = mu, phi = phi)
  a[is.na(a)] <- -99
  b[is.na(b)] <- -99
  u <- stats::runif(n = length(y), min = a, max = b)
  u[u == -99] <- NA
  stats::qnorm(u)
}

qres_truncated_nbinom1 <- function(object, y, mu, ...) {
  theta <- get_pars(object)
  phi <- exp(theta[["ln_phi"]])
  if (is_delta(object)) phi <- phi[2]
  a <- ptruncated_nbinom1(y - 1, mu = mu, phi = phi)
  b <- ptruncated_nbinom1(y, mu = mu, phi = phi)
  a[is.na(a)] <- -99
  b[is.na(b)] <- -99
  u <- stats::runif(n = length(y), min = a, max = b)
  u[u == -99] <- NA
  stats::qnorm(u)
}

qres_pois <- function(object, y, mu, ...) {
  a <- stats::ppois(y - 1, mu)
  b <- stats::ppois(y, mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}

is_delta <- function(object) {
  isTRUE(object$family$delta)
}

qres_gamma <- function(object, y, mu, ...) {
  theta <- get_pars(object)
  phi <- exp(theta[["ln_phi"]])
  if (is_delta(object)) phi <- phi[2]
  s1 <- phi
  s2 <- mu / s1
  u <- stats::pgamma(q = y, shape = s1, scale = s2)
  stats::qnorm(u)
}

qres_gamma_mix <- function(object, y, mu, ...) {
  cli_abort("Randomized quantile residuals for this family are not implemented yet")
  # theta <- get_pars(object)
  # p_extreme <- plogis(theta[["logit_p_extreme"]])
  # phi <- exp(theta[["ln_phi"]])
  # if (is_delta(object)) phi <- phi[2]
  # ratio <- exp(theta[["log_ratio_mix"]])
  # s1 <- phi
  # s2 <- mu / s1
  # s3 <- (ratio * mu) / s1
  # u <- stats::pgamma(q = y, shape = s1, scale = (1-p_extreme)*s2 + p_extreme*s3) # this looks wrong
  # stats::qnorm(u)
}

qres_nbinom2_mix <- function(object, y, mu, ...) {
  cli_abort("Randomized quantile residuals for this family are not implemented yet")
  theta <- get_pars(object)
  p_extreme <- plogis(theta[["logit_p_extreme"]])
  phi <- exp(theta[["ln_phi"]])
  if (is_delta(object)) phi <- phi[2]
  ratio <- exp(theta[["log_ratio_mix"]])
  a <- stats::pnbinom(y - 1, size = phi, mu = (1-p_extreme)*mu + p_extreme*ratio*mu)
  b <- stats::pnbinom(y, size = phi, mu = (1-p_extreme)*mu + p_extreme*ratio*mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}

qres_lognormal_mix <- function(object, y, mu, ...) {
  cli_abort("Randomized quantile residuals for this family are not implemented yet")
  theta <- get_pars(object)
  p_extreme <- plogis(theta[["logit_p_extreme"]])
  dispersion <- exp(theta[["ln_phi"]])
  if (is_delta(object)) dispersion <- dispersion[2]
  ratio <- exp(theta[["log_ratio_mix"]])
  u <- stats::plnorm(q = y, meanlog = log((1-p_extreme)*mu + p_extreme*ratio*mu) - (dispersion^2) / 2, sdlog = dispersion)
  stats::qnorm(u)
}

qres_gaussian <- function(object, y, mu, ...) {
  theta <- get_pars(object)
  dispersion <- exp(theta[["ln_phi"]])
  u <- stats::pnorm(q = y, mean = mu, sd = dispersion)
  stats::qnorm(u)
}

qres_lognormal <- function(object, y, mu, ...) {
  theta <- get_pars(object)
  dispersion <- exp(theta[["ln_phi"]])
  if (is_delta(object)) dispersion <- dispersion[2]
  u <- stats::plnorm(q = y, meanlog = log(mu) - (dispersion^2) / 2, sdlog = dispersion)
  stats::qnorm(u)
}

# https://en.wikipedia.org/wiki/Location%E2%80%93scale_family
pt_ls <- function(q, df, mu, sigma) stats::pt((q - mu) / sigma, df)

qres_student <- function(object, y, mu, ...) {
  theta <- get_pars(object)
  dispersion <- exp(theta[["ln_phi"]])
  df <- exp(theta[["ln_student_df"]]) + 1
  u <- pt_ls(q = y, df = df, mu = mu, sigma = dispersion)
  stats::qnorm(u)
}

qres_beta <- function(object, y, mu, ...) {
  theta <- get_pars(object)
  phi <- exp(theta[["ln_phi"]])
  s1 <- mu * phi
  s2 <- (1 - mu) * phi
  u <- stats::pbeta(q = y, shape1 = s1, shape2 = s2)
  stats::qnorm(u)
}

# Modified from https://github.com/chjackson/flexsurv/blob/c765344fa798868036b841481fce2ea4d009d85e/src/gengamma.h#L63
pgengamma <- function(q, mean, sigma, .Q, lower.tail = TRUE, log.p = FALSE) {
  if (.Q != 0) {
    y <- log(q)
    beta <- .Q / sigma
    k <- .Q^-2
    mu <- log(mean) - lgamma((k * beta + 1) / beta) + lgamma(k) + log(k) / beta # Need to convert from mean to 'location' mu
    w <- (y - mu) / sigma
    expnu <- exp(.Q * w) * k

    if (.Q > 0) {
      stats::pgamma(q = expnu, shape = k, rate = 1, lower.tail = lower.tail, log.p = log.p)
    } else {
      stats::pgamma(q = expnu, shape = k, rate = 1, lower.tail = !lower.tail, log.p = log.p)
    }
  } else {
    # use lnorm with correction for Q = 0
    stats::plnorm(q = q, meanlog = log(mean) - ((sigma^2) / 2), sdlog = sigma)
  }
}

qres_gengamma <- function(object, y, mu, ...) {
  theta <- get_pars(object)
  .Q <- theta$gengamma_Q
  sigma <- exp(theta$ln_phi)
  if (is_delta(object)) sigma <- sigma[2]
  u <- pgengamma(q = y, mean = mu, sigma = sigma, .Q = .Q)
  stats::qnorm(u)
}

#' Residuals method for sdmTMB models
#'
#' See the residual-checking vignette: `browseVignettes("sdmTMB")` or [on the
#' documentation
#' site](https://sdmTMB.github.io/sdmTMB/articles/residual-checking.html).
#' See notes about types of residuals in 'Details' section below.
#'
#' @param object An [sdmTMB()] model.
#' @param type Residual type. See details.
#' @param model Which delta/hurdle model component?
#' @param mcmc_samples A vector of MCMC samples of the linear predictor in link
#'   space. See the `predict_mle_mcmc()` function in the
#'   \href{https://github.com/sdmTMB/sdmTMBextra}{sdmTMBextra} package.
#' @param qres_func A custom quantile residuals function. Function should take
#'   the arguments `object, y, mu, ...` and return a vector of length
#'   `length(y)`.
#' @param ... Passed to custom `qres_func` function. Unused.
#' @export
#' @importFrom stats predict
#' @details
#'
#' **Randomized quantile residuals:**
#'
#' `mle-mvn`, `mle-eb`, and `mle-mcmc` are all implementations of
#' randomized quantile residuals (Dunn & Smyth 1996), which are also known as
#' probability integral transform (PIT) residuals (Smith 1985). If the data are
#' consistent with model assumptions, these residuals should be distributed as
#' normal(0, 1). Randomization is added to account for integer or binary
#' response observations. For example, for a Poisson observation likelihood with
#' observations `y` and mean predictions `mu`, we would create randomized
#' quantile residuals as:
#'
#' ```
#' a <- ppois(y - 1, mu)
#' b <- ppois(y, mu)
#' u <- runif(n = length(y), min = a, max = b)
#' qnorm(u)
#' ```
#'
#' **Types of residuals:**
#'
#' Acronyms:
#' - EB: Empirical Bayes
#' - MCMC: Markov chain Monte Carlo
#' - MLE: Maximum Likelihood Estimate
#' - MVN: Multivariate normal
#'
#' **`mle-mvn`**: Fixed effects are held at their MLEs and random effects are
#' taken from a single approximate posterior sample. The "approximate" part
#' refers to the sample being taken from the random effects' assumed MVN
#' distribution. In practice, the sample is obtained based on the mode and
#' Hessian of the random effects taking advantage of sparsity in the Hessian for
#' computational efficiency. This sample is taken with `obj$MC()`, where `obj`
#' is the \pkg{TMB} object created with `TMB::MakeADFun()`. See Waagepetersen
#' (2006) and the description in the source code for the internal \pkg{TMB}
#' function `TMB:::oneSamplePosterior()`. Residuals are converted to randomized
#' quantile residuals as described above.
#'
#' **`mle-eb`**: Fixed effects are held at their MLEs and random effects are
#' taken as their EB estimates. These used to be the default residuals in
#' \pkg{sdmTMB} (and were called `mle-laplace`). They are available for
#' backwards compatibility and for research purposes but they are *not*
#' recommended for checking goodness of fit. Residuals are converted to
#' randomized quantile residuals as described above.
#'
#' **`mle-mcmc`**: Fixed effects are held at their MLEs and random effects are
#' taken from a single posterior sample obtained with MCMC. These are an
#' excellent option since they make no assumption about the distribution of the
#' random effects (compared to the `mle-mvn` option) but can be slow to obtain.
#' See Waagepetersen (2006) and Thygesen et al. (2017). Residuals are converted
#' to randomized quantile residuals as described above.
#'
#' See the \href{https://github.com/sdmTMB/sdmTMBextra}{\pkg{sdmTMBextra}}
#' package for the function `predict_mle_mcmc()`, which can generate the MCMC
#' samples to pass to the `mcmc_samples` argument. Ideally MCMC is run until
#' convergence and then the last iteration can be used for residuals.
#' The defaults may not be sufficient for many models.
#'
#' **`response`**: These are simple observed minus predicted residuals.
#'
#' **`pearson`**: These are Pearson residuals: response residuals scaled by the
#' standard deviation. If weights are present, the residuals are then
#' multiplied by sqrt(weights).
#'
#' @references
#' Dunn, P.K. & Smyth, G.K. (1996). Randomized Quantile Residuals. Journal of
#' Computational and Graphical Statistics, 5, 236–244.
#'
#' Smith, J.Q. (1985). Diagnostic checks of non-standard time series models.
#' Journal of Forecasting, 4, 283–291.
#'
#' Waagepetersen, R. (2006). A simulation-based goodness-of-fit test for random
#' effects in generalized linear mixed models. Scandinavian Journal of
#' Statistics, 33(4), 721-731.
#'
#' Thygesen, U.H., Albertsen, C.M., Berg, C.W., Kristensen, K., and Nielsen, A.
#' 2017. Validation of ecological state space models using the Laplace
#' approximation. Environ Ecol Stat 24(2): 317–339.
#' \doi{10.1007/s10651-017-0372-4}
#'
#' Rufener, M.-C., Kristensen, K., Nielsen, J.R., and Bastardie, F. 2021.
#' Bridging the gap between commercial fisheries and survey data to model the
#' spatiotemporal dynamics of marine species. Ecological Applications. e02453.
#' \doi{10.1002/eap.2453}
#'
#' @return A vector of residuals. Note that randomization from any single
#' random effect posterior sample and from any randomized quantile routines
#' will result in different residuals with each call. It is suggested to **set
#' a randomization seed** and to not go "fishing" for the perfect residuals or
#' to present all inspected residuals.
#'
#' @seealso [simulate.sdmTMB()], [dharma_residuals()]
#' @examples
#'
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 10)
#' fit <- sdmTMB(
#'   present ~ as.factor(year) + poly(depth, 2),
#'   data = pcod_2011, mesh = mesh,
#'   family = binomial()
#' )
#'
#' # the default "mle-mvn" residuals use fixed effects at their MLE and a
#' # single sample from the approximate random effect posterior:
#' set.seed(9283)
#' r <- residuals(fit, type = "mle-mvn")
#' qqnorm(r)
#' abline(0, 1)
#'
#' # response residuals will be not be normally distributed unless
#' # the family is Gaussian:
#' r <- residuals(fit, type = "response")
#' qqnorm(r)
#' abline(0, 1)
#'
#' # "mle-eb" are quick but are not expected to be N(0, 1); not recommended:
#' set.seed(2321)
#' r <- residuals(fit, type = "mle-eb")
#' qqnorm(r)
#' abline(0, 1)
#'
#' # see also "mle-mcmc" residuals with the help of the sdmTMBextra package
#' # we can fake them here by taking a single sample from the joint precision
#' # matrix and pretending they are MCMC samples:
#' set.seed(82728)
#' p <- predict(fit, nsim = 1) # pretend these are from sdmTMBextra::predict_mle_mcmc()
#' r <- residuals(fit, mcmc_samples = p)
#' qqnorm(r)
#' abline(0, 1)
residuals.sdmTMB <- function(object,
                             type = c("mle-mvn", "mle-eb", "mle-mcmc", "response", "pearson", "deviance"),
                             model = c(1, 2),
                             mcmc_samples = NULL,
                             qres_func = NULL,
                             ...) {
  type_was_missing <- missing(type)
  type <- match.arg(type[[1]], choices = c("mle-mvn", "mle-laplace", "mle-eb", "mle-mcmc", "response", "pearson", "deviance"))

  # retrieve function that called this:
  sys_calls <- unlist(lapply(sys.calls(), deparse))
  visreg_call <- any(grepl("setupV", substr(sys_calls, 1, 7)))
  if (type == "mle-laplace") type <- "mle-eb"
  model_missing <- FALSE
  if (identical(model, c(1, 2))) model_missing <- TRUE
  model <- as.integer(model[[1]])
  if ("visreg_model" %in% names(object)) {
    model <- object$visreg_model
  }
  # need to re-attach environment if in fresh session
  reinitialize(object)

  fam <- object$family$family
  nd <- NULL
  est_column <- "est"
  linkinv <- object$family$linkinv
  if (isTRUE(object$family$delta)) {
    fam <- fam[[model]]
    linkinv <-  object$family[[model]]$linkinv
    nd <- object$data
    est_column <- if (model == 1L) "est1" else "est2"
  }
  if(fam %in% c("truncated_nbinom1", "truncated_nbinom2")){
    linkinv <- function(eta){exp(eta)}
  } # for residuals, use untruncated mean
  if (is.null(qres_func)) {
    res_func <- switch(fam,
      gaussian = qres_gaussian,
      binomial = qres_binomial,
      betabinomial = qres_betabinomial,
      tweedie  = qres_tweedie,
      Beta     = qres_beta,
      Gamma    = qres_gamma,
      nbinom2  = qres_nbinom2,
      nbinom1  = qres_nbinom1,
      truncated_nbinom2  = qres_truncated_nbinom2,
      truncated_nbinom1  = qres_truncated_nbinom1,
      poisson  = qres_pois,
      student  = qres_student,
      lognormal  = qres_lognormal,
      gamma_mix = qres_gamma_mix,
      lognormal_mix = qres_lognormal_mix,
      nbinom2_mix = qres_nbinom2_mix,
      gengamma = qres_gengamma,
      cli_abort(paste(fam, "not yet supported."))
    )
  } else {
    res_func <- qres_func
  }

  if (!"offset" %in% names(object)) cli_abort("This model appears to have been fit with an older sdmTMB.")
  if (type %in% c("mle-eb", "response", "pearson")) {
    mu <- linkinv(predict(object, newdata = object$data, offset = object$offset)[[est_column]]) # not newdata = NULL
  } else if (type == "mvn-laplace") {
    mu <- linkinv(predict(object, nsim = 1L, model = model, offset = object$offset)[, 1L, drop = TRUE])
  } else if (type == "mle-mcmc") {
    if (is.null(mcmc_samples)) {
      msg <- c("As of sdmTMB 0.3.0, `mcmc_samples` must be supplied to use `type = 'mle-mcmc'`.",
        "See ?sdmTMBextra::predict_mle_mcmc after installing",
        "remotes::install_github('sdmTMB/sdmTMBextra')")
      cli_abort(msg)
    }
    mcmc_samples <- as.numeric(mcmc_samples)
    assert_that(length(mcmc_samples) == nrow(object$data))
    mu <- linkinv(mcmc_samples)
  } else if (type == "mle-mvn") {
    ## see TMB:::oneSamplePosterior()

    if (is.null(object$tmb_random)) {
      params <- object$tmb_obj$env$last.par.best
    } else {
      params <- .one_sample_posterior(object)
    }
    pred <- predict(
      object,
      newdata = object$data,
      mcmc_samples = matrix(params, ncol = 1L),
      model = model[[1L]],
      nsim = 1L,
      offset = object$offset
    )
    mu <- linkinv(pred[, 1L, drop = TRUE])
  } else if (type == "deviance") {
    # if (is_delta(object)) {
    #   cli_abort("Deviance residuals not implemented for delta models")
    # }
    if (!is_delta(object)) {
      if (object$family$family=="binomial" && !all(object$tmb_data$size == 1)) {
        cli_abort("Deviance residuals not implemented for binomial family with size > 1.")
      }
    }
    report <- object$tmb_obj$report(object$tmb_obj$env$last.par.best)
    resids <- report$devresid[, model, drop = TRUE]
    if (all(resids == 0)) {
      cli_abort("Deviance residuals not implemented for this family.")
    }
    return(resids)
  } else {
    cli_abort("residual type not implemented")
  }

  y <- object$response
  y <- y[, model, drop = TRUE] # in case delta
  # e.g., visreg, prediction has already removed NA mu:
  if (sum(is.na(y)) > 0 && length(mu) < length(y)) {
    y <- y[!is.na(y)]
  }

  # for binomial proportion with weights = N:
  size <- object$tmb_data$size
  prop_binomial <- !all(size == 1)

  if (type == "response") {
    if (!prop_binomial) r <- y - mu else r <- y / size - mu
  } else if (type == "mle-eb" || type == "mle-mvn") {
    r <- res_func(object, y, mu, .n = size, ...)
  } else if (type == "mle-mcmc") {
    r <- res_func(object, y, mu, .n = size, ...)
  } else if (type == "pearson") {
    if (is.null(v <- family(object)$variance)) {
      cli_abort(c("Variance function undefined for family;",
        "cannot compute Pearson residuals"))
    }
    # FIXME: add sigma function, for now just binomial
    if (length(formals(v)) > 1)
      cli_abort("sdmTMB currently only supports variance functions with 1 argument")
    # vv <- switch(length(formals(v)),
    #   v(fitted(object)),
    #   v(fitted(object), sigma(object)),
    #   stop("variance function should take 1 or 2 arguments"))
    vv <- v(mu)
    wts <- if (prop_binomial) size else object$tmb_data$weights_i
    if (!prop_binomial) r <- y - mu else r <- y / size - mu
    r <- r / sqrt(vv)
    if (!is.null(wts)) r <- r * sqrt(wts)
  } else {
    cli_abort("residual type not implemented")
  }
  if (isTRUE(object$family$delta) && is.null(mcmc_samples) && model_missing) {
    cli_inform(paste0("These are residuals for delta model component ", model,
      ". Use the `model` argument to select the other component."))
  }
  as.vector(r)
}

# from:
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor
check_overdisp <- function(object) {
  rdf <- stats::df.residual(object)
  rp <- stats::residuals(object, type = "pearson")
  pearson_chisq <- sum(rp^2)
  prat <- pearson_chisq / rdf
  pval <- stats::pchisq(pearson_chisq, df = rdf, lower.tail = FALSE)
  data.frame(chisq = pearson_chisq, ratio = prat, rdf = rdf, p = pval)
}

# return full set of parameter vector with the
# random effects sampled from the implied MVN posterior and the
# fixed effects at their MLEs
.one_sample_posterior <- function(object) {
  tmp <- object$tmb_obj$env$MC(n = 1L, keep = TRUE, antithetic = FALSE)
  re_samp <- as.vector(attr(tmp, "samples"))
  lp <- object$tmb_obj$env$last.par.best
  p <- numeric(length(lp))
  fe <- object$tmb_obj$env$lfixed()
  re <- object$tmb_obj$env$lrandom()
  p[re] <- re_samp
  p[fe] <- lp[fe]
  p
}
