# Reviewer, lane wt-priorform: does frm(prior = validate_prior(pl, ...))
# fit identically to frm(prior = pl)? Also: which parameters get which
# density, read off resolve_prior_input(), in both spellings.
#   Rscript dev/priorform-rev-roundtrip.R      (lane build; seed 20260916)
.libPaths(c("C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(20260916)
ng <- 15; per <- 12; n <- ng * per
d <- data.frame(g = factor(rep(1:ng, each = per)), x = runif(n, 0, 3),
                z = rnorm(n), time = factor(rep(1:4, length.out = n)))
u <- rnorm(ng, 0, .5)
d$y <- 1 + 0.5 * d$x + u[d$g] + (0.3 * rnorm(ng))[d$g] * d$z + rnorm(n)
d$ynl <- 3 * exp(-0.6 * d$x) + u[d$g] * 0.3 + rnorm(n, 0, 0.2)
d$y1 <- d$y + rnorm(n); d$y2 <- 0.5 * d$z + u[d$g] + rnorm(n)
lat <- 0.8 * d$x - 0.6 * d$z + u[d$g] + rlogis(n)
d$yo <- cut(lat, c(-Inf, 0, 1.2, 2.4, Inf), labels = FALSE)
d$yb <- pmin(pmax(plogis(0.3 * d$x - 0.5 + u[d$g] * 0.5 + rnorm(n, 0, .3)),
                  0.01), 0.99)

cases <- list(
  beta_sd_classwide = list(bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(),
    set_prior("normal(0, 0.1)", class = "sd") + set_prior("normal(0, 2)", class = "b")),
  beta_sd_group_then_class = list(bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(),
    set_prior("normal(0, 5)", class = "sd", group = "g") +
      set_prior("normal(0, 0.05)", class = "sd")),
  beta_sd_class_then_group = list(bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(),
    set_prior("normal(0, 0.05)", class = "sd") +
      set_prior("normal(0, 5)", class = "sd", group = "g")),
  beta_sd_phi_only = list(bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(),
    set_prior("normal(0, 0.05)", class = "sd", group = "g", dpar = "phi")),
  b_coef_then_class = list(bf(y ~ x + z) + gaussian(),
    set_prior("normal(3, 0.01)", coef = "x") + set_prior("normal(0, 0.01)", class = "b")),
  b_class_then_coef = list(bf(y ~ x + z) + gaussian(),
    set_prior("normal(0, 0.01)", class = "b") + set_prior("normal(3, 0.01)", coef = "x")),
  bounds_then_density = list(bf(y ~ x + z) + gaussian(),
    set_prior("", coef = "x", lb = 0.9) + set_prior("normal(0, 1)", coef = "x")),
  density_then_bounds = list(bf(y ~ x + z) + gaussian(),
    set_prior("normal(0, 1)", coef = "x") + set_prior("", coef = "x", lb = 0.9)),
  class_bound_then_coef_density = list(bf(y ~ x + z) + gaussian(),
    set_prior("", class = "b", lb = 0.9) + set_prior("normal(0, 1)", coef = "z")),
  mv = list(mvbf(bf(y1 ~ x + (1 | p | g)), bf(y2 ~ z + (1 | p | g))),
    set_prior("normal(0, 0.1)", class = "b", resp = "y1") +
      set_prior("normal(0, 0.1)", class = "sd", resp = "y2") +
      set_prior("lkj(4)", class = "cor")),
  nl = list(bf(ynl ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE) + gaussian(),
    set_prior("normal(2, 0.1)", nlpar = "a") + set_prior("", nlpar = "b", lb = 0, ub = 5) +
      set_prior("normal(0, 0.1)", class = "sd", nlpar = "a")),
  smooth = list(bf(y ~ s(x, k = 5) + (1 | g)) + gaussian(),
    set_prior("normal(0, 0.5)", class = "b") + set_prior("normal(0, 0.2)", class = "sd")),
  gp_theta = list(bf(y ~ gp(x)) + gaussian(),
    set_prior("normal(0, 0.3)", class = "theta")),
  ordinal = list(bf(yo ~ x + z + (1 | g)) + cumulative(),
    set_prior("normal(0, 1)", class = "Intercept") + set_prior("normal(0, 0.5)", class = "b", coef = "z")),
  sratio_cs = list(bf(yo ~ x + cs(z)) + sratio(),
    set_prior("normal(0, 0.3)", class = "b")),
  cor_slope = list(bf(y ~ x + (z | g)) + gaussian(),
    set_prior("lkj(10)", class = "cor") + set_prior("normal(0, 0.2)", class = "sd", group = "g")),
  sigma_natural = list(bf(y ~ x, sigma ~ z) + gaussian(),
    set_prior("normal(0, 0.1)", class = "b", dpar = "sigma") +
      set_prior("normal(0, 0.5)", class = "Intercept", dpar = "sigma"))
)
fitsum <- function(f) list(ll = as.numeric(logLik(f)), par = f$opt$par,
                           obj = f$opt$objective)
entries <- function(f, pl) {
  e <- frmtmb:::resolve_prior_input(f, pl)
  paste(sort(c(vapply(e$entries, function(z) paste0(z$comp, paste(z$idx, collapse = ","), "=",
                                             paste(unlist(z$dist), collapse = "/")), ""),
               paste0("lb:", names(e$lower), "=", e$lower),
               paste0("ub:", names(e$upper), "=", e$upper))), collapse = " ")
}
for (nm in names(cases)) {
  cs <- cases[[nm]]
  res <- tryCatch({
    vp <- validate_prior(cs[[2]], cs[[1]], data = d)
    f1 <- suppressWarnings(suppressMessages(frm(cs[[1]], data = d, prior = cs[[2]])))
    f2 <- suppressWarnings(suppressMessages(frm(cs[[1]], data = d, prior = vp)))
    same <- identical(fitsum(f1), fitsum(f2))
    e1 <- entries(f1, cs[[2]]); e2 <- entries(f1, vp)
    sprintf("fit identical=%s  entries identical=%s  dll=%.3g%s", same,
            identical(e1, e2),
            as.numeric(logLik(f1)) - as.numeric(logLik(f2)),
            if (!identical(e1, e2)) paste0("\n    direct : ", e1, "\n    roundtr: ", e2) else "")
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(nm, ":", res, "\n")
}
# the table itself for the order case, as a user would read it
print(validate_prior(cases$b_coef_then_class[[2]], cases$b_coef_then_class[[1]], data = d))
print(validate_prior(cases$beta_sd_group_then_class[[2]],
                     cases$beta_sd_group_then_class[[1]], data = d))
