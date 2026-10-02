# Build a model matrix and a compact specification of its encoding.
#
# The specification holds the evaluated model-frame terms (whose `predvars`
# retain fitted transformation bases such as `poly()`), factor levels,
# contrasts, and the expected columns. `.apply_design()` reproduces the same
# encoding on new rows without refitting to the original data. `name` labels
# the design in error messages.
.make_design <- function(formula, data, name) {
  mf <- stats::model.frame(stats::delete.response(stats::terms(formula)),
    data, na.action = stats::na.pass)
  terms <- attr(mf, "terms")
  X <- stats::model.matrix(terms, mf)
  if (anyNA(X)) {
    cli_abort("Missing values were produced in the {name} design matrix.")
  }
  spec <- list(
    name = name,
    terms = terms,
    xlevels = stats::.getXlevels(terms, mf),
    contrasts = attr(X, "contrasts"),
    columns = colnames(X)
  )
  list(X = X, spec = spec)
}

# Apply a specification from `.make_design()` to new rows. Returns a matrix
# with one row per row of `newdata` and the specification's columns in order.
# A caller may shrink `spec$columns` to a subset of the full matrix columns.
.apply_design <- function(spec, newdata) {
  mf <- tryCatch(
    stats::model.frame(spec$terms, newdata, xlev = spec$xlevels,
      na.action = stats::na.pass),
    error = function(e) {
      cli_abort(c(
        "Could not build the {spec$name} design matrix for `newdata`.",
        "x" = "{conditionMessage(e)}",
        "i" = "Check that factor levels and columns match those used in fitting."
      ), call = NULL)
    }
  )
  X <- stats::model.matrix(spec$terms, mf, contrasts.arg = spec$contrasts)
  if (nrow(X) != nrow(newdata)) {
    cli_abort("Internal error: the {spec$name} design matrix has {nrow(X)} row{?s} but `newdata` has {nrow(newdata)}.")
  }
  missing_cols <- setdiff(spec$columns, colnames(X))
  if (length(missing_cols)) {
    cli_abort(c(
      "The {spec$name} design matrix for `newdata` does not match the fitted model.",
      "x" = "Missing column{?s}: {missing_cols}."
    ))
  }
  X <- X[, spec$columns, drop = FALSE]
  if (anyNA(X)) {
    cli_abort(c(
      "Missing values are not allowed in variables used by {spec$name} in `newdata`.",
      "x" = "Rows with missing values: {which(!stats::complete.cases(X))}."
    ))
  }
  X
}
