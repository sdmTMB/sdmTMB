# Multi-family post-fit capability policy.
#
# This is the single source of truth for which methods support models fit with
# a named list of families (`n_f > 1`). Supported methods are adapted for
# row-wise families and tested. Every other method must call
# `.check_family_capability()` immediately and abort, rather than partially
# interpreting multi-family metadata or silently operating on the first family.
.multi_family_supported_methods <- c(
  # inspection
  "print", "tidy", "coef", "confint", "fixef", "ranef", "vcov",
  "logLik", "AIC", "nobs", "df.residual", "family",
  # prediction and simulation
  "predict", "fitted", "simulate", "dharma_residuals",
  # derived quantities (prediction rows must select a single target family
  # where a combined index is required; see the multi-family vignette)
  "get_index", "get_index_split", "get_cog", "get_eao", "get_weighted_average",
  # field plots
  "plot_anisotropy", "plot_anisotropy2", "calculate_anisotropy_components"
)

# Abort early if `object` is a multi-family fit and `method` is not supported.
# Returns the family_spec (invisibly) so supported paths can reuse it without a
# second lookup. `caller` is the method name as it should appear in the error
# message, e.g. "`residuals()`". `info` optionally adds an "i" hint bullet.
.check_family_capability <- function(object, method, caller = NULL, info = NULL) {
  if (is.null(caller)) {
    caller <- paste0("`", method, "()`")
  }
  family_spec <- .object_family_spec(object, caller = caller)
  if (.family_spec_is_multi_family(family_spec) &&
      !method %in% .multi_family_supported_methods) {
    msg <- "{caller} is not yet supported for multi-family models."
    if (!is.null(info)) {
      msg <- c(msg, "i" = info)
    }
    cli_abort(msg)
  }
  invisible(family_spec)
}
