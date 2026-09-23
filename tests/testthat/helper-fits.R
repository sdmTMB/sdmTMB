# Every 4th cell of `qcs_grid` in each direction: 452 rows instead of 7,314.
qcs_grid_small <- qcs_grid[qcs_grid$X %% 8 == 0 & qcs_grid$Y %% 8 == 0, ]

# Memoise an expensive fit so test blocks in a file can share it. The fit runs
# on first use, so it still sits behind each block's skips.
fit_once <- function(f) {
  val <- NULL
  function() {
    if (is.null(val)) val <<- f()
    val
  }
}

# Full-size simulation-recovery and Monte Carlo tests.
skip_slow <- function() {
  skip_if_not(nzchar(Sys.getenv("SDMTMB_SLOW_TESTS")), "slow test")
}
