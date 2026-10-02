# Reference-fit case table. Sourced by run.R after sdmTMB is loaded, with
# `rf_dir` set to the reference-fits directory.
#
# Each case is a list with:
#   args      arguments to sdmTMB() (backend and getsd are set by run.R)
#   newdata   prediction data
#   pred_offset  newdata column used as the prediction offset, if any
#   extras    optional extra metrics: "re_form_na", "index", "cog"
#   tol       optional per-metric tolerance overrides, e.g. list(nll = c(rel = 1e-5))
#   allow_sanity  sanity() checks allowed to fail when recording

read_fixture <- function(file) {
  utils::read.csv(file.path(rf_dir, "data", file), stringsAsFactors = FALSE)
}

fam <- read_fixture("family.csv")
fam_nd <- read_fixture("family-newdata.csv")

sp <- read_fixture("spatial.csv")
sp$g <- factor(sp$g)
sp_nd <- read_fixture("spatial-newdata.csv")
sp_nd$g <- factor(sp_nd$g, levels = levels(sp$g))

# Nonlocal covariate on a dense grid covering every mesh vertex and year.
# x_surface() is the same deterministic surface used in make-data.R.
x_surface <- function(X, Y, year) {
  sin(X / 2.5) + cos(Y / 3) * 0.8 + 0.25 * (year - 2013) * cos(X / 4)
}
nl_grid <- expand.grid(X = seq(-1, 11, 0.5), Y = seq(-1, 11, 0.5), year = 2011:2015)
nl_grid$x <- x_surface(nl_grid$X, nl_grid$Y, nl_grid$year)

mesh <- make_mesh(sp, c("X", "Y"), mesh = fmesher::fm_rcdt_2d_inla(
  loc = as.matrix(read_fixture("mesh-loc.csv")),
  tv = as.matrix(read_fixture("mesh-tv.csv"))
))

ohio <- as.data.frame(ohio_df)
ohio$log_pop <- log(ohio$pop)
ohio_edges <- read_fixture("ohio-edges.csv")
ohio_domain <- make_areal_domain(
  igraph::graph_from_data_frame(ohio_edges, directed = FALSE,
    vertices = data.frame(name = sort(unique(ohio$county)))),
  space_column = "county"
)

cases <- list()
add_case <- function(name, ..., newdata, pred_offset = NULL, extras = character(),
                     tol = list(), allow_sanity = character()) {
  if (!is.null(cases[[name]])) stop("Duplicate case name: ", name)
  cases[[name]] <<- list(args = list(...), newdata = newdata,
    pred_offset = pred_offset, extras = extras, tol = tol,
    allow_sanity = allow_sanity)
}

# Families and links (non-spatial) --------------------------------------------

family_links <- list(
  gaussian = c("identity", "log", "inverse"),
  student = c("identity", "log", "inverse"),
  Gamma = c("inverse", "identity", "log"),
  lognormal = c("identity", "log", "inverse"),
  gengamma = c("identity", "log", "inverse"),
  poisson = c("log", "identity"),
  nbinom2 = "log",
  nbinom1 = "log",
  truncated_nbinom2 = "log",
  truncated_nbinom1 = "log",
  tweedie = c("log", "identity"),
  gamma_mix = c("log", "identity", "inverse"),
  lognormal_mix = c("log", "identity", "inverse"),
  nbinom2_mix = "log",
  binomial = c("logit", "cloglog"),
  Beta = "logit",
  ordbeta = "logit",
  censored_poisson = "log"
)

make_family <- function(family, link) {
  if (family == "student") return(suppressMessages(student(link = link, df = NULL)))
  get(family, mode = "function")(link = link)
}
# Identity and inverse links need starting values near the generating
# intercept (b = 0 gives an invalid mean), as a user would supply.
start_b <- function(link) {
  switch(link, identity = c(2, 0), inverse = c(0.5, 0), NULL)
}
for (f in names(family_links)) {
  for (l in family_links[[f]]) {
    y <- paste0("y_", f, "_", l)
    start <- if (!is.null(start_b(l))) list(b_j = start_b(l)) else list()
    control <- if (f == "censored_poisson") {
      sdmTMBcontrol(censored_upper = fam$upr_censored_poisson)
    } else {
      sdmTMBcontrol(start = start)
    }
    add_case(paste0("family/", f, "_", l),
      formula = stats::as.formula(paste(y, "~ x")), data = fam,
      family = make_family(f, l), spatial = "off", control = control,
      newdata = fam_nd)
  }
}
for (l in c("logit", "cloglog")) {
  add_case(paste0("family/betabinomial_", l),
    formula = stats::as.formula(paste0("cbind(y_betabinomial_", l,
      ", size - y_betabinomial_", l, ") ~ x")),
    data = fam, family = betabinomial(link = l), spatial = "off",
    newdata = fam_nd)
}
add_case("family/binomial_trials_cbind",
  formula = cbind(y_binomial_trials, size - y_binomial_trials) ~ x,
  data = fam, family = binomial(), spatial = "off", newdata = fam_nd)
add_case("family/binomial_trials_weights",
  formula = I(y_binomial_trials / size) ~ x, weights = fam$size,
  data = fam, family = binomial(), spatial = "off", newdata = fam_nd)

offset_links <- list(gaussian = "identity", poisson = "log", Gamma = "inverse",
  binomial = c("logit", "cloglog"))
for (f in names(offset_links)) {
  for (l in offset_links[[f]]) {
    add_case(paste0("offset/", f, "_", l),
      formula = stats::as.formula(paste0("y_", f, "_", l, "_off ~ x")),
      data = fam, offset = "off", family = make_family(f, l), spatial = "off",
      newdata = fam_nd, pred_offset = "off")
  }
}

# Delta families (non-spatial) ---------------------------------------------------

delta_families <- list(
  Gamma = delta_gamma(),
  lognormal = delta_lognormal(),
  gengamma = delta_gengamma(),
  truncated_nbinom2 = delta_truncated_nbinom2(),
  truncated_nbinom1 = delta_truncated_nbinom1(),
  gamma_mix = delta_gamma_mix(),
  lognormal_mix = delta_lognormal_mix(),
  Beta = delta_beta(),
  Gamma_cloglog = delta_gamma(link1 = "cloglog"),
  Gamma_inverse = delta_gamma(link2 = "inverse"),
  Gamma_identity = delta_gamma(link2 = "identity"),
  lognormal_cloglog = delta_lognormal(link1 = "cloglog")
)
delta_start <- list(
  Gamma_inverse = list(b_j2 = c(0.4, 0)),
  Gamma_identity = list(b_j2 = c(2.5, 0))
)
for (nm in names(delta_families)) {
  start <- if (is.null(delta_start[[nm]])) list() else delta_start[[nm]]
  add_case(paste0("delta/", nm),
    formula = stats::as.formula(paste0("y_delta_", nm, " ~ x")), data = fam,
    family = delta_families[[nm]], spatial = "off", newdata = fam_nd,
    control = sdmTMBcontrol(start = start))
}
delta_offset <- list(
  Gamma = delta_gamma(),
  lognormal = delta_lognormal(),
  truncated_nbinom2 = delta_truncated_nbinom2()
)
for (nm in names(delta_offset)) {
  add_case(paste0("delta/", nm, "_offset"),
    formula = stats::as.formula(paste0("y_delta_", nm, "_off ~ x")), data = fam,
    offset = "off", family = delta_offset[[nm]], spatial = "off",
    newdata = fam_nd, pred_offset = "off")
}
delta_pl <- list(
  Gamma = delta_gamma(type = "poisson-link"),
  lognormal = delta_lognormal(type = "poisson-link"),
  gengamma = delta_gengamma(type = "poisson-link"),
  lognormal_mix = delta_lognormal_mix(type = "poisson-link")
)
for (nm in names(delta_pl)) {
  add_case(paste0("delta/poisson_link_", nm),
    formula = stats::as.formula(paste0("y_deltapl_", nm, " ~ x")), data = fam,
    family = delta_pl[[nm]], spatial = "off", newdata = fam_nd)
}
for (nm in c("Gamma", "lognormal")) {
  add_case(paste0("delta/poisson_link_", nm, "_offset"),
    formula = stats::as.formula(paste0("y_deltapl_", nm, "_off ~ x")), data = fam,
    offset = "off", family = delta_pl[[nm]], spatial = "off",
    newdata = fam_nd, pred_offset = "off")
}

# Random-field configurations ---------------------------------------------------

sp_case <- function(name, formula, ..., data = sp, newdata = sp_nd) {
  add_case(paste0("fields/", name), formula = formula, data = data,
    mesh = if (identical(data, sp)) mesh else make_mesh(data, c("X", "Y"), mesh = mesh$mesh),
    ..., newdata = newdata)
}
f_year <- y_gauss ~ 0 + as.factor(year) + x
sp_case("spatial_only", y_gauss ~ x, extras = "re_form_na")
sp_case("st_iid_no_spatial", f_year, time = "year", spatial = "off",
  spatiotemporal = "iid")
sp_case("st_iid", f_year, time = "year", spatiotemporal = "iid",
  extras = c("re_form_na", "index", "cog"))
sp_case("st_ar1", y_gauss ~ x, time = "year", spatiotemporal = "ar1")
sp_case("st_ar1_no_spatial", y_gauss ~ x, time = "year", spatial = "off",
  spatiotemporal = "ar1")
sp_case("st_rw", y_gauss ~ x, time = "year", spatiotemporal = "rw")
sp_case("st_iid_reml", f_year, time = "year", spatiotemporal = "iid", reml = TRUE)
sp_case("share_range_false", y_range ~ 0 + as.factor(year) + x, time = "year",
  spatiotemporal = "iid", share_range = FALSE)
sp_case("anisotropy", y_aniso ~ x, anisotropy = TRUE)
sp_case("spatial_varying", y_svc ~ x, spatial_varying = ~ 0 + x,
  extras = "re_form_na")
for (type in c("rw", "rw0", "ar1")) {
  sp_case(paste0("time_varying_", type), if (type == "rw") y_tv ~ 0 else y_tv ~ 1 + z,
    time = "year", spatiotemporal = "off", time_varying = ~ 1 + z,
    time_varying_type = type)
}
sp_case("extra_time", y_gauss ~ x, time = "year", spatiotemporal = "ar1",
  extra_time = 2013L, data = sp[sp$year != 2013L, ])
sp_case("iid_intercept", y_re ~ x + (1 | g), extras = "re_form_na")
sp_case("iid_correlated_slope", y_re ~ x + (1 + x | g))
sp_case("smoother", y_smooth ~ s(z))
# The breakpoint objective has a kink, so gradients at the optimum can be large.
sp_case("breakpt", y_breakpt ~ breakpt(z), allow_sanity = c("nlminb_ok", "gradients_ok"))
sp_case("logistic", y_logistic ~ logistic(z))
sp_case("dispformula", y_gauss ~ x, dispformula = ~ z)
sp_case("priors", y_gauss ~ x, priors = sdmTMBpriors(
  matern_s = pc_matern(range_gt = 1, sigma_lt = 2),
  b = normal(c(0, 0), c(5, 5))
))
sp_case("bayesian", y_gauss ~ x, bayesian = TRUE, priors = sdmTMBpriors(
  matern_s = pc_matern(range_gt = 1, sigma_lt = 2)
))
sp_case("nb2_st_iid_offset", y_nb2 ~ x, time = "year", spatiotemporal = "iid",
  family = nbinom2(), offset = "off", pred_offset = "off",
  extras = c("index"))
sp_case("nb2_st_ar1_offset", y_nb2 ~ x, time = "year", spatiotemporal = "ar1",
  family = nbinom2(), offset = "off", pred_offset = "off")
sp_case("gamma_spatial", y_gamma ~ x, family = Gamma(link = "log"))
sp_case("multi_family", y_mf ~ x,
  family = list(gauss = gaussian(), pois = poisson()),
  distribution_column = "dist")
nonlocal_terms <- list(diffusion = ~ diffusion(x), time_lag = ~ time_lag(x),
  combined = ~ diffusion(x) + time_lag(x))
for (nm in names(nonlocal_terms)) {
  sp_case(paste0("nonlocal_", nm), y_nl ~ x, time = "year", spatial = "off",
    spatiotemporal = "off", nonlocal_formula = nonlocal_terms[[nm]],
    nonlocal_data = nl_grid)
}

# Delta models with spatial fields: shared and separate structure ------------

d_case <- function(name, ..., family = delta_gamma()) {
  sp_case(paste0("delta_", name), ..., family = family, offset = "off",
    pred_offset = "off")
}
d_case("both_iid", y_delta ~ x + z, time = "year", spatiotemporal = "iid",
  extras = c("index", "re_form_na"))
d_case("formula_list", list(y_delta ~ x, y_delta ~ z), time = "year",
  spatiotemporal = "iid")
d_case("spatial_list", list(y_delta ~ x, y_delta ~ z), time = "year",
  spatial = list("on", "off"), spatiotemporal = list("iid", "off"))
d_case("st_list_ar1", list(y_delta ~ x, y_delta ~ z), time = "year",
  spatiotemporal = list("off", "ar1"))
d_case("share_range_list", list(y_delta ~ x, y_delta ~ z), time = "year",
  spatiotemporal = "iid", share_range = list(FALSE, TRUE))
d_case("anisotropy", list(y_delta_aniso ~ x, y_delta_aniso ~ z), anisotropy = TRUE)
d_case("spatial_varying", list(y_delta_svc ~ x, y_delta_svc ~ z),
  spatial_varying = ~ 0 + x)
d_case("time_varying", list(y_delta_tv ~ 0 + x, y_delta_tv ~ 0 + z), time = "year",
  spatiotemporal = "off", time_varying = ~ 1)
d_case("iid_intercepts", list(y_delta_re ~ x + (1 | g), y_delta_re ~ z + (1 | g)))
d_case("dispformula", list(y_delta ~ x, y_delta ~ z), dispformula = ~ z)
d_case("poisson_link", list(y_deltapl ~ x, y_deltapl ~ z), time = "year",
  spatiotemporal = list("iid", "off"), family = delta_lognormal(type = "poisson-link"),
  extras = "index")

# Areal SAR/CAR ---------------------------------------------------------------------

for (model in c("sar", "car")) {
  add_case(paste0("areal/", model, "_spatial"),
    formula = cases ~ pct_male, data = ohio, mesh = ohio_domain,
    spatial_model = model, family = poisson(), offset = "log_pop",
    newdata = ohio, pred_offset = "log_pop")
}
# Ohio has too little spatiotemporal variation for SAR (sigma_E -> 0), so only
# CAR gets a spatiotemporal case.
for (model in "car") {
  add_case(paste0("areal/", model, "_st_iid"),
    formula = cases ~ 0 + as.factor(year) + pct_male, data = ohio,
    mesh = ohio_domain, spatial_model = model, time = "year",
    spatiotemporal = "iid", family = poisson(), offset = "log_pop",
    newdata = ohio, pred_offset = "log_pop")
}
