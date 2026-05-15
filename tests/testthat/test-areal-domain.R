make_areal_test_graph <- function(W) {
  igraph::graph_from_adjacency_matrix(
    as.matrix(W),
    mode = "directed",
    weighted = TRUE,
    diag = TRUE
  )
}

test_that("make_areal_domain normalizes igraph input", {
  W <- Matrix::Matrix(
    c(
      0, 2, 0,
      1, 0, 1,
      0, 3, 0
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(W) <- colnames(W) <- c("a", "b", "c")

  d <- make_areal_domain(make_areal_test_graph(W), space_column = "region")

  expect_s3_class(d, "sdmTMBareal")
  expect_identical(d$n_s, 3L)
  expect_identical(d$unit_names, c("a", "b", "c"))
  expect_equal(as.matrix(d$W_raw), as.matrix(W))
  expect_equal(as.numeric(Matrix::rowSums(d$W)), c(1, 1, 1))
  expect_equal(unname(Matrix::diag(d$W)), c(0, 0, 0))
})

test_that("make_areal_domain handles islands", {
  W <- Matrix::Matrix(
    c(
      0, 1, 0,
      1, 0, 0,
      0, 0, 0
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(W) <- colnames(W) <- c("a", "b", "c")

  expect_warning(
    d <- make_areal_domain(make_areal_test_graph(W), space_column = "region"),
    "islands"
  )
  expect_equal(as.numeric(Matrix::rowSums(d$W)), c(1, 1, 0))
})

test_that("make_areal_domain works with named igraph input", {
  g <- igraph::make_empty_graph(n = 3, directed = FALSE)
  igraph::V(g)$name <- c("a", "b", "c")
  g <- igraph::add_edges(g, c("a", "b", "b", "c"))

  d <- make_areal_domain(g, space_column = "region")

  expect_s3_class(d, "sdmTMBareal")
  expect_identical(d$unit_names, c("a", "b", "c"))
  expect_equal(as.numeric(Matrix::rowSums(d$W)), c(1, 1, 1))
})

test_that("make_areal_domain works with labelled sf polygon input", {
  skip_if_not_installed("sf")

  grid <- sf::st_make_grid(
    sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 1))),
    n = c(2, 1),
    square = TRUE
  )
  poly <- sf::st_sf(area = c("left", "right"), geometry = grid)
  d <- make_areal_domain(poly, id_column = "area")

  expect_s3_class(d, "sdmTMBareal")
  expect_identical(d$unit_names, c("left", "right"))
  expect_identical(d$space_column, "area")
  expect_equal(dim(d$W), c(2L, 2L))
  expect_equal(as.numeric(Matrix::rowSums(d$W)), c(1, 1))
})

test_that("make_areal_grid overlays point data on an existing sf grid", {
  skip_if_not_installed("sf")

  grid <- sf::st_make_grid(
    sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 2))),
    n = c(2, 2),
    square = TRUE
  )
  dat <- data.frame(
    x = c(0.25, 1.25, 0.25, 1.25),
    y = c(0.25, 0.25, 1.25, 1.25),
    z = 1:4
  )

  out <- make_areal_grid(
    data = dat,
    xy_cols = c("x", "y"),
    spatial_domain = grid,
    space_column = "cell"
  )

  expect_true("cell" %in% names(out$data))
  expect_s3_class(out$domain, "sdmTMBareal")
  expect_s3_class(out$grid, "sf")
  expect_equal(out$domain$n_s, 4L)
  expect_equal(sort(unique(out$data$cell)), sort(out$domain$unit_names))
  expect_equal(as.numeric(Matrix::rowSums(out$domain$W)), rep(1, 4))
})

test_that("make_areal_grid can create and clip a grid from an sf boundary", {
  skip_if_not_installed("sf")

  boundary <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 2))))
  dat <- data.frame(
    x = c(0.25, 1.25, 0.25, 1.25),
    y = c(0.25, 0.25, 1.25, 1.25)
  )

  out <- make_areal_grid(
    data = dat,
    xy_cols = c("x", "y"),
    spatial_domain = boundary,
    n = c(2, 2),
    square = TRUE
  )

  expect_s3_class(out$domain, "sdmTMBareal")
  expect_equal(out$domain$n_s, 4L)
  expect_equal(nrow(out$grid), 4L)
  expect_true("grid_cell" %in% names(out$data))
})

test_that("prepare_spatial_domain validates data memberships", {
  W <- Matrix::Matrix(
    c(
      0, 1,
      1, 0
    ),
    nrow = 2,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(W) <- colnames(W) <- c("a", "b")
  d <- make_areal_domain(make_areal_test_graph(W), space_column = "region")

  expect_error(
    prepare_spatial_domain(
      mesh = d,
      data = data.frame(region = c("a", "x")),
      mesh_missing = FALSE,
      anisotropy = FALSE,
      covariate_diffusion = NULL,
      spatial_model = "sar"
    ),
    "no match in the areal domain"
  )

  expect_identical(d$unit_names, c("a", "b"))
})

test_that("prepare_spatial_domain validates required membership column and missing values", {
  W <- Matrix::Matrix(
    c(
      0, 1,
      1, 0
    ),
    nrow = 2,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(W) <- colnames(W) <- c("a", "b")
  d <- make_areal_domain(make_areal_test_graph(W), space_column = "region")

  expect_error(
    prepare_spatial_domain(
      mesh = d,
      data = data.frame(x = 1:2),
      mesh_missing = FALSE,
      anisotropy = FALSE,
      covariate_diffusion = NULL,
      spatial_model = "sar"
    ),
    "areal domain was not found in `data`"
  )
  expect_error(
    prepare_spatial_domain(
      mesh = d,
      data = data.frame(region = c("a", NA)),
      mesh_missing = FALSE,
      anisotropy = FALSE,
      covariate_diffusion = NULL,
      spatial_model = "sar"
    ),
    "contains missing values"
  )
})

test_that("make_areal_domain validates malformed graph inputs", {

  g_unnamed <- igraph::make_ring(3, directed = FALSE)
  expect_error(
    make_areal_domain(g_unnamed, space_column = "region"),
    "vertex names"
  )

  expect_error(
    make_areal_domain(matrix(c(0, 1, 1, 0), nrow = 2), space_column = "region"),
    "named igraph object or an sf/sfc polygon object"
  )

  g_duplicate <- igraph::make_empty_graph(n = 2, directed = FALSE)
  igraph::V(g_duplicate)$name <- c("a", "a")
  expect_error(
    make_areal_domain(g_duplicate, space_column = "region"),
    "must be unique"
  )

  W_diag <- Matrix::Matrix(c(1, 1, 1, 0), nrow = 2, byrow = TRUE, sparse = TRUE)
  rownames(W_diag) <- colnames(W_diag) <- c("a", "b")
  expect_error(
    make_areal_domain(make_areal_test_graph(W_diag), space_column = "region"),
    "self-neighbor loops"
  )
})

test_that("areal projection helpers map memberships", {
  W <- Matrix::Matrix(
    c(
      0, 1,
      1, 0
    ),
    nrow = 2,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(W) <- colnames(W) <- c("a", "b")
  d <- make_areal_domain(make_areal_test_graph(W), space_column = "region")

  idx <- domain_obs_index(d, data.frame(region = c("b", "a", "b")))
  expect_identical(idx, c(2L, 1L, 2L))

  A <- areal_projection_matrix(d, data.frame(region = c("b", "a", "b")))
  expect_equal(dim(A), c(3, 2))
  expect_equal(as.vector(as.matrix(A)), as.vector(matrix(c(
    0, 1,
    1, 0,
    0, 1
  ), nrow = 3, byrow = TRUE)))
})

test_that("prepare_spatial_domain returns areal domain pieces", {
  W <- Matrix::Matrix(
    c(
      0, 1,
      1, 0
    ),
    nrow = 2,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(W) <- colnames(W) <- c("a", "b")
  dat <- data.frame(region = c("a", "b", "a"))
  d <- make_areal_domain(make_areal_test_graph(W), space_column = "region")

  out <- prepare_spatial_domain(
    mesh = d,
    data = dat,
    mesh_missing = FALSE,
    anisotropy = FALSE,
    covariate_diffusion = NULL,
    spatial_model = "sar"
  )

  expect_identical(out$type, "areal")
  expect_identical(out$n_s, 2L)
  expect_equal(dim(out$A_st), c(3L, 2L))
  expect_equal(out$A_spatial_index, 0:2)
  expect_s4_class(out$W_ss, "dgCMatrix")
  expect_equal(dim(out$W_ss), c(2L, 2L))
  expect_equal(as.matrix(out$W_ss), as.matrix(d$W))
  expect_identical(out$normalize_in_r, 0L)
  expect_identical(out$barrier, 0L)
  expect_identical(out$anisotropy, 0L)
})

test_that("prepare_spatial_domain returns raw adjacency for CAR", {
  W <- Matrix::Matrix(
    c(
      0, 1, 0,
      1, 0, 1,
      0, 1, 0
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(W) <- colnames(W) <- c("a", "b", "c")
  dat <- data.frame(region = c("a", "b", "c", "a"))
  d <- make_areal_domain(make_areal_test_graph(W), space_column = "region")

  out <- prepare_spatial_domain(
    mesh = d,
    data = dat,
    mesh_missing = FALSE,
    anisotropy = FALSE,
    covariate_diffusion = NULL,
    spatial_model = "car"
  )

  expect_identical(out$type, "areal")
  expect_identical(out$spatial_model, "car")
  expect_equal(as.matrix(out$W_ss), as.matrix(d$W_raw))
})

test_that("prepare_spatial_domain can return raw adjacency for SAR", {
  W <- Matrix::Matrix(
    c(
      0, 1, 0,
      1, 0, 1,
      0, 1, 0
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(W) <- colnames(W) <- c("a", "b", "c")
  dat <- data.frame(region = c("a", "b", "c", "a"))
  d <- make_areal_domain(make_areal_test_graph(W), space_column = "region")

  out <- prepare_spatial_domain(
    mesh = d,
    data = dat,
    mesh_missing = FALSE,
    anisotropy = FALSE,
    covariate_diffusion = NULL,
    spatial_model = "sar",
    sar_weight_style = "raw"
  )

  expect_identical(out$spatial_model, "sar")
  expect_equal(as.matrix(out$W_ss), as.matrix(d$W_raw))
})

test_that("prepare_spatial_domain validates unsupported areal options", {
  W <- Matrix::Matrix(
    c(
      0, 1,
      1, 0
    ),
    nrow = 2,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(W) <- colnames(W) <- c("a", "b")
  dat <- data.frame(region = c("a", "b", "a"), x = c(1, 2, 3))
  d <- make_areal_domain(make_areal_test_graph(W), space_column = "region")

  expect_error(
    prepare_spatial_domain(
      mesh = d,
      data = dat,
      mesh_missing = FALSE,
      share_range = FALSE,
      anisotropy = FALSE,
      covariate_diffusion = NULL,
      spatial_model = "sar"
    ),
    "share_range = FALSE"
  )

  expect_error(
    prepare_spatial_domain(
      mesh = d,
      data = dat,
      mesh_missing = FALSE,
      anisotropy = TRUE,
      covariate_diffusion = NULL,
      spatial_model = "sar"
    ),
    "anisotropy.*not supported"
  )

  expect_error(
    prepare_spatial_domain(
      mesh = d,
      data = dat,
      mesh_missing = FALSE,
      anisotropy = FALSE,
      covariate_diffusion = list(dummy = 1),
      spatial_model = "sar"
    ),
    "covariate_diffusion.*not supported"
  )

  expect_error(
    prepare_spatial_domain(
      mesh = d,
      data = dat["x"],
      mesh_missing = FALSE,
      anisotropy = FALSE,
      covariate_diffusion = NULL,
      spatial_model = "sar"
    ),
    "areal domain was not found in `data`"
  )

  pri <- sdmTMBpriors(matern_s = pc_matern(range_gt = 2, sigma_lt = 2))
  expect_error(
    prepare_spatial_domain(
      mesh = d,
      data = dat,
      mesh_missing = FALSE,
      anisotropy = FALSE,
      priors = pri,
      covariate_diffusion = NULL,
      spatial_model = "sar"
    ),
    "PC Matern priors.*not supported"
  )

  d_barrier <- d
  d_barrier$spde_barrier <- TRUE
  expect_error(
    prepare_spatial_domain(
      mesh = d_barrier,
      data = dat,
      mesh_missing = FALSE,
      anisotropy = FALSE,
      covariate_diffusion = NULL,
      spatial_model = "sar"
    ),
    "Barrier models are not supported"
  )

  expect_error(
    prepare_spatial_domain(
      mesh = d,
      data = dat,
      mesh_missing = FALSE,
      anisotropy = FALSE,
      covariate_diffusion = NULL,
      experimental = list(epsilon_model = ~x),
      spatial_model = "sar"
    ),
    "epsilon_model.*not supported"
  )
})

test_that("set_limits applies default bounds for SAR rho", {
  tmb_obj <- list(par = c(logit_rho_sar = 0))
  lim <- set_limits(
    tmb_obj = tmb_obj,
    lower = list(),
    upper = list()
  )
  expect_equal(
    unname(lim$lower["logit_rho_sar"]),
    stats::qlogis((-0.999 + 1) / 2)
  )
  expect_equal(
    unname(lim$upper["logit_rho_sar"]),
    stats::qlogis((0.999 + 1) / 2)
  )
})

test_that("set_limits applies default bounds for CAR alpha", {
  tmb_obj <- list(par = c(logit_rho_sar = 0))
  lim <- set_limits(
    tmb_obj = tmb_obj,
    lower = list(),
    upper = list(),
    spatial_model = 2L
  )
  expect_equal(unname(lim$lower["logit_rho_sar"]), -Inf)
  expect_equal(unname(lim$upper["logit_rho_sar"]), stats::qlogis(0.999))
})

build_predict_stub_with_areal_domain <- function(getsd = FALSE) {
  skip_if_not_installed("TMB")

  set.seed(11)
  fit_dat <- data.frame(
    y = stats::rnorm(24),
    x = stats::rnorm(24),
    X = seq(0, 23) / 10,
    Y = rep(c(0, 1, 2), each = 8)
  )
  mesh <- make_mesh(fit_dat, xy_cols = c("X", "Y"), cutoff = 0.4)
  fit <- sdmTMB(
    y ~ x,
    data = fit_dat,
    mesh = mesh,
    family = gaussian(),
    control = sdmTMBcontrol(getsd = getsd),
    silent = TRUE
  )

  areal_data <- data.frame(region = rep(c("a", "b", "c"), each = 8L))
  W <- Matrix::Matrix(
    c(
      0, 1, 0,
      1, 0, 1,
      0, 1, 0
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(W) <- colnames(W) <- c("a", "b", "c")
  domain <- make_areal_domain(make_areal_test_graph(W), space_column = "region")

  fit$data$region <- areal_data$region
  fit$spde <- domain
  fit$tmb_data$spatial_model <- 1L
  fit$tmb_data$no_spatial <- 0L
  list(fit = fit, domain = domain)
}

test_that("areal predict uses one-hot projection matrix for spatial predictions", {
  obj <- build_predict_stub_with_areal_domain()
  fit <- obj$fit
  domain <- obj$domain
  nd <- data.frame(
    x = c(-0.25, 0.25, 0.5, -0.5),
    region = c("a", "c", "b", "a")
  )
  td <- predict(fit, newdata = nd, return_tmb_data = TRUE)

  expect_s4_class(td$proj_mesh, "dgCMatrix")
  expect_equal(dim(td$proj_mesh), c(nrow(nd), domain$n_s))
  expect_equal(as.numeric(Matrix::rowSums(td$proj_mesh)), rep(1, nrow(nd)))
  expect_equal(td$proj_spatial_index, seq_len(nrow(nd)) - 1L)
})

test_that("areal predict allows population predictions without space column", {
  obj <- build_predict_stub_with_areal_domain()
  fit <- obj$fit
  domain <- obj$domain
  nd <- data.frame(x = seq(-1, 1, length.out = 4))
  td <- predict(fit, newdata = nd, re_form = NA, return_tmb_data = TRUE)

  expect_equal(dim(td$proj_mesh), c(nrow(nd), domain$n_s))
  expect_identical(Matrix::nnzero(td$proj_mesh), 0L)
  expect_equal(td$proj_spatial_index, seq_len(nrow(nd)) - 1L)
})

test_that("areal predict errors on unknown areal units", {
  obj <- build_predict_stub_with_areal_domain()
  fit <- obj$fit

  expect_error(
    predict(fit, newdata = data.frame(x = 0, region = "z"), return_tmb_data = TRUE),
    "not found in the domain"
  )
})

test_that("areal tidy/print/sanity reporting paths avoid Matérn assumptions", {
  obj <- suppressWarnings(build_predict_stub_with_areal_domain(getsd = TRUE))
  fit <- obj$fit

  td <- suppressWarnings(tidy(fit, effects = "ran_pars", silent = TRUE))
  expect_false("range" %in% td$term)
  expect_true("sigma_O" %in% td$term)

  txt <- suppressWarnings(paste(capture.output(print(fit)), collapse = "\n"))
  expect_match(txt, "SAR field scale")
  expect_false(grepl("Mat..rn range", txt))

  s <- suppressMessages(suppressWarnings(sanity(fit)))
  expect_true(is.list(s))
  expect_true("range_ok" %in% names(s))
})

test_that("spread_sims uses areal field-scale transform when areal flag is set", {
  skip_if_not_installed("TMB")

  set.seed(21)
  dat <- data.frame(
    y = stats::rnorm(60),
    x = stats::rnorm(60),
    X = stats::runif(60),
    Y = stats::runif(60)
  )
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.2)
  fit_spde <- sdmTMB(
    y ~ x,
    data = dat,
    mesh = mesh,
    family = gaussian(),
    control = sdmTMBcontrol(get_joint_precision = TRUE),
    silent = TRUE
  )

  fit_areal <- fit_spde
  fit_areal$tmb_data$spatial_model <- 1L
  fit_areal$spde <- structure(list(), class = c("sdmTMBareal", "sdmTMBdomain"))

  set.seed(42)
  sims_spde <- spread_sims(fit_spde, nsim = 8)
  set.seed(42)
  sims_areal <- spread_sims(fit_areal, nsim = 8)

  ln_kappa <- log(sqrt(8) / sims_spde$range)
  ln_tau_O <- -log(sims_spde$sigma_O * sqrt(4 * pi)) - ln_kappa
  expected_sigma_O_areal <- exp(-ln_tau_O)

  expect_false("range" %in% names(sims_areal))
  expect_equal(sims_areal$sigma_O, expected_sigma_O_areal, tolerance = 1e-10)
})

build_areal_smoke_domain <- function() {
  dat <- data.frame(
    region = rep(c("a", "b", "c"), each = 4L),
    y = c(-1.2, -0.8, -1.0, -0.6, -0.1, 0.2, 0.0, 0.3, 0.7, 1.0, 0.8, 1.2)
  )
  W <- Matrix::Matrix(
    c(
      0, 1, 0,
      1, 0, 1,
      0, 1, 0
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(W) <- colnames(W) <- c("a", "b", "c")
  list(
    data = dat,
    domain = make_areal_domain(make_areal_test_graph(W), space_column = "region")
  )
}

test_that("minimal areal do_fit FALSE model builds TMB object", {
  skip_if_not_installed("TMB")

  smoke <- build_areal_smoke_domain()
  fit <- sdmTMB(
    y ~ 1,
    data = smoke$data,
    mesh = smoke$domain,
    spatial_model = "sar",
    family = gaussian(),
    do_fit = FALSE,
    silent = TRUE
  )

  expect_s3_class(fit, "sdmTMB")
  expect_true("logit_rho_sar" %in% names(fit$tmb_obj$par))
})

test_that("areal domains require explicit SAR or CAR spatial_model", {
  skip_if_not_installed("TMB")

  smoke <- build_areal_smoke_domain()
  expect_error(
    sdmTMB(
      y ~ 1,
      data = smoke$data,
      mesh = smoke$domain,
      family = gaussian(),
      do_fit = FALSE,
      silent = TRUE
    ),
    "spatial_model = \"spde\""
  )

  dat <- data.frame(y = stats::rnorm(8), X = seq_len(8), Y = rep(1:2, each = 4))
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.1)
  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial_model = "sar",
      family = gaussian(),
      do_fit = FALSE,
      silent = TRUE
    ),
    "requires an areal domain"
  )
})

test_that("minimal Gaussian spatial-only areal SAR model fits", {
  skip_if_not_installed("TMB")

  smoke <- build_areal_smoke_domain()
  fit <- sdmTMB(
    y ~ 1,
    data = smoke$data,
    mesh = smoke$domain,
    spatial_model = "sar",
    family = gaussian(),
    control = sdmTMBcontrol(getsd = FALSE, newton_loops = 0),
    silent = TRUE
  )

  expect_s3_class(fit, "sdmTMB")
  expect_true("logit_rho_sar" %in% names(fit$tmb_obj$par))
})

test_that("tidy reports SAR rho on transformed scale for areal models", {
  skip_if_not_installed("TMB")

  smoke <- build_areal_smoke_domain()
  fit <- sdmTMB(
    y ~ 1,
    data = smoke$data,
    mesh = smoke$domain,
    spatial_model = "sar",
    family = gaussian(),
    control = sdmTMBcontrol(newton_loops = 0),
    silent = TRUE
  )
  td <- tidy(fit, "ran_pars", silent = TRUE)
  rho <- td[td$term == "rho_sar", , drop = FALSE]

  expect_false("range" %in% td$term)
  expect_equal(nrow(rho), 1L)
  expect_true(all(rho$estimate >= -1 & rho$estimate <= 1))
  expect_true(all(rho$conf.low >= -1 & rho$conf.high <= 1))
})

test_that("areal SAR can use raw weights", {
  smoke <- build_areal_smoke_domain()
  fit <- sdmTMB(
    y ~ 1,
    data = smoke$data,
    mesh = smoke$domain,
    spatial_model = "sar",
    family = gaussian(),
    control = sdmTMBcontrol(
      getsd = FALSE,
      newton_loops = 0,
      sar_weight_style = "raw"
    ),
    do_fit = FALSE,
    silent = TRUE
  )

  expect_equal(as.matrix(fit$tmb_data$W_ss), as.matrix(smoke$domain$W_raw))
})

test_that("minimal Gaussian spatial-only areal CAR model fits", {
  skip_if_not_installed("TMB")

  smoke <- build_areal_smoke_domain()
  fit <- sdmTMB(
    y ~ 1,
    data = smoke$data,
    mesh = smoke$domain,
    spatial_model = "car",
    family = gaussian(),
    control = sdmTMBcontrol(getsd = FALSE, newton_loops = 0),
    silent = TRUE
  )

  expect_s3_class(fit, "sdmTMB")
  expect_true("logit_rho_sar" %in% names(fit$tmb_obj$par))
  expect_identical(fit$tmb_data$spatial_model, 2L)
})

test_that("tidy reports CAR alpha on transformed scale for areal models", {
  skip_if_not_installed("TMB")

  smoke <- build_areal_smoke_domain()
  fit <- sdmTMB(
    y ~ 1,
    data = smoke$data,
    mesh = smoke$domain,
    spatial_model = "car",
    family = gaussian(),
    control = sdmTMBcontrol(newton_loops = 0),
    silent = TRUE
  )
  td <- tidy(fit, "ran_pars", silent = TRUE)
  alpha <- td[td$term == "alpha_car", , drop = FALSE]

  expect_false("range" %in% td$term)
  expect_false("rho_sar" %in% td$term)
  expect_equal(nrow(alpha), 1L)
  expect_true(all(alpha$estimate >= 0 & alpha$estimate <= 1))
  expect_true(all(alpha$conf.low >= 0 & alpha$conf.high <= 1))
})

test_that("areal share_range FALSE errors before TMB object construction", {
  skip_if_not_installed("TMB")

  smoke <- build_areal_smoke_domain()
  expect_error(
    sdmTMB(
      y ~ 1,
      data = smoke$data,
      mesh = smoke$domain,
      spatial_model = "sar",
      family = gaussian(),
      share_range = FALSE,
      spatial = "off",
      do_fit = FALSE,
      silent = TRUE
    ),
    "share_range = FALSE"
  )
})

test_that("no-spatial areal models do not estimate SAR rho", {
  skip_if_not_installed("TMB")

  smoke <- build_areal_smoke_domain()
  fit <- sdmTMB(
    y ~ 1,
    data = smoke$data,
    mesh = smoke$domain,
    spatial_model = "sar",
    family = gaussian(),
    spatial = "off",
    do_fit = FALSE,
    silent = TRUE
  )

  expect_s3_class(fit, "sdmTMB")
  expect_false("logit_rho_sar" %in% names(fit$tmb_obj$par))
})
