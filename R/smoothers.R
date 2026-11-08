# from brms:::rm_wsp()
rm_wsp <- function (x) {
  out <- gsub("[ \t\r\n]+", "", x, perl = TRUE)
  dim(out) <- dim(x)
  out
}
# from brms:::all_terms()
all_terms <- function (x) {
  if (!length(x)) {
    return(character(0))
  }
  if (!inherits(x, "terms")) {
    x <- terms(stats::as.formula(x))
  }
  rm_wsp(attr(x, "term.labels"))
}

get_smooth_terms <- function(terms) {
  x1 <- grep("s\\(", terms)
  x2 <- grep("t2\\(", terms)
  c(x1, x2)
}

parse_smoothers <- function(formula, data, knots = NULL, newdata = NULL, basis_prev = NULL) {
  terms <- all_terms(formula)
  smooth_i <- get_smooth_terms(terms)
  basis <- list()
  basis_out <- list()
  Zs <- list()
  Xs <- list()
  zs_basis <- integer(0) # basis index for each penalized block in `Zs`
  labels <- list()
  classes <- list()
  if (length(smooth_i) > 0) {
    has_smooths <- TRUE
    smterms <- terms[smooth_i]
    ns <- 0
    ns_Xf <- 0
    for (i in seq_along(smterms)) {
      if (grepl('bs\\=\\"re', smterms[i])) cli_abort("bs = 're' is not currently supported for smooths")
      if (grepl('fx\\=T', smterms[i])) cli_abort("fx = TRUE is not currently supported for smooths")
      # if (grepl('m\\=3', smterms[i])) cli_abort("m > 2 is not currently supported for smooths")
      obj <- eval(str2expression(smterms[i]))
      labels[[i]] <- obj$label
      classes[[i]] <- attr(obj, "class")
      if (is.null(newdata)) {
        basis_out[[i]] <- basis[[i]] <- mgcv::smoothCon(
          object = obj, data = data,
          knots = knots, absorb.cons = TRUE, modCon = 3,
          diagonal.penalty = FALSE
        )
      } else {
        basis[[i]] <- basis_prev[[i]] # predicting on new data
      }
      for (j in seq_along(basis[[i]])) { # elements > 1 with `by` terms
        ns_Xf <- ns_Xf + 1
        rasm <- mgcv::smooth2random(basis[[i]][[j]], names(data), type = 2)
        if (!is.null(newdata)) {
          rasm <- s2rPred(basis[[i]][[j]], rasm, newdata)
        }
        for (k in seq_along(rasm$rand)) { # elements > 1 with if s(x, y) or t2()
          ns <- ns + 1
          Zs[[ns]] <- rasm$rand[[k]]
          zs_basis[ns] <- ns_Xf
        }
        Xs[[ns_Xf]] <- rasm$Xf
      }
    }
    sm_dims <- unlist(lapply(Zs, ncol))
    xf_dims <- vapply(Xs, ncol, integer(1L)) # unpenalized (null space) columns per basis
    Xs <- do.call(cbind, Xs) # combine 'em all into one design matrix
    b_smooth_start <- c(0, cumsum(sm_dims)[-length(sm_dims)])
  } else {
    has_smooths <- FALSE
    sm_dims <- 0L
    xf_dims <- integer(0)
    b_smooth_start <- 0L
    Xs <- matrix(nrow = 0L, ncol = 0L)
  }
  list(Xs = Xs, Zs = Zs, has_smooths = has_smooths, labels = labels,
    classes = classes, basis_out = basis_out,
    sm_dims = sm_dims, b_smooth_start = b_smooth_start,
    xf_dims = xf_dims, zs_basis = zs_basis)
}

# from mgcv docs ?mgcv::smooth2random
s2rPred <- function(sm, re, data) {
  ## Function to aid prediction from smooths represented as type==2
  ## random effects. re must be the result of smooth2random(sm,...,type=2).
  if (!all(sm$term %in% colnames(data))) {
    cli_abort(paste("A smoother term is missing from 'newdata':",
      sm$term[!sm$term %in% colnames(data)]))
  }
  X <- mgcv::PredictMat(sm, data) ## get prediction matrix for new data
  ## transform to r.e. parameterization
  if (!is.null(re$trans.U)) {
    X <- X %*% re$trans.U
  }
  X <- t(t(X) * re$trans.D)
  ## re-order columns according to random effect re-ordering...
  X[, re$rind] <- X[, re$pen.ind != 0]
  ## re-order penalization index in same way
  pen.ind <- re$pen.ind
  pen.ind[re$rind] <- pen.ind[pen.ind > 0]
  ## start return object...
  r <- list(rand = list(), Xf = X[, which(re$pen.ind == 0), drop = FALSE])
  for (i in seq_along(re$rand)) { ## loop over random effect matrices
    r$rand[[i]] <- X[, which(pen.ind == i), drop = FALSE]
    attr(r$rand[[i]], "s.label") <- attr(re$rand[[i]], "s.label")
  }
  names(r$rand) <- names(re$rand)
  r
}
## use function to obtain prediction random and fixed effect matrices
## for first 10 elements of 'dat'. Then confirm that these match the
## first 10 rows of the original model matrices, as they should...
# r <- s2rPred(sm,re,dat[1:10,])
# range(r$Xf-re$Xf[1:10,])
# range(r$rand[[1]]-re$rand[[1]][1:10,])
