.is_multi_family_model <- function(object) {
  !is.null(object$tmb_data) &&
    !is.null(object$tmb_data$multi_family) &&
    isTRUE(object$tmb_data$multi_family == 1L)
}

.has_multi_family_delta <- function(object) {
  .is_multi_family_model(object) &&
    !is.null(object$tmb_data$delta_family_e) &&
    any(object$tmb_data$delta_family_e == 1L)
}

.is_delta_like_model <- function(object) {
  if (.is_multi_family_model(object)) {
    .has_multi_family_delta(object)
  } else {
    isTRUE(object$family$delta)
  }
}

.is_regular_delta_model <- function(object) {
  !.is_multi_family_model(object) &&
    isTRUE(object$family$delta)
}
