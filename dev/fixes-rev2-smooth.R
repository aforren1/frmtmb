# Reviewer of lane fixes, re-check: the s() basis switch and the
# smooth_fx_units() optimizer scaling, base against lane, over gamSim
# designs and smooth shapes the lane did not run. One record per fit.
#   Rscript dev/fixes-rev2-smooth.R <lib> <seed_from> <seed_to> <out.rds>
a <- commandArgs(TRUE)
lib <- a[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
seeds <- seq(as.integer(a[2]), as.integer(a[3]))
gs <- function(eg, seed, n = 200, re = TRUE, ...) {
  set.seed(seed)
  d <- suppressMessages(mgcv::gamSim(eg = eg, n = n, verbose = FALSE, ...))
  if (eg == 2) d <- d$data else d$z <- runif(nrow(d))
  d$g <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
  d$g8 <- factor(sample(letters[1:8], nrow(d), TRUE))
  d$gr <- factor(rep(1:20, length.out = nrow(d)))
  if (re) d$y <- d$y + 0.6 * stats::rnorm(20)[d$gr]
  d
}
cases <- list(
  F1 = list(6, bf(y ~ s(x1) + s(x2)), gaussian(), y ~ s(x1) + s(x2)),
  F2 = list(6, bf(y ~ s(x1, bs = "cr", k = 6)), gaussian(),
            y ~ s(x1, bs = "cr", k = 6)),
  F3 = list(6, bf(y ~ s(x1, by = g) + g), gaussian(), y ~ s(x1, by = g) + g),
  F4 = list(6, bf(y ~ s(x1, by = z)), gaussian(), y ~ s(x1, by = z)),
  F5 = list(6, bf(y ~ s(x1, x2)), gaussian(), y ~ s(x1, x2)),
  F6 = list(6, bf(y ~ t2(x1, x2)), gaussian(), y ~ t2(x1, x2)),
  F7 = list(6, bf(y ~ s(x0) + s(x1) + s(x2) + s(x3)), gaussian(),
            y ~ s(x0) + s(x1) + s(x2) + s(x3)),
  F8 = list(6, bf(y ~ x1 + s(gr, bs = "re")), gaussian(),
            y ~ x1 + s(gr, bs = "re")),
  F9 = list(6, bf(y ~ s(x1, g, bs = "fs", k = 5)), gaussian(),
            y ~ s(x1, g, bs = "fs", k = 5)),
  F10 = list(6, bf(y ~ g8 + s(x1, by = g8, k = 5)), gaussian(),
             y ~ g8 + s(x1, by = g8, k = 5)),
  F11 = list(6, bf(y ~ s(x1), sigma ~ s(x2)), gaussian(), NULL),
  F12 = list(6, bf(y ~ s(x1) + s(x2) + (1 | gr)), gaussian(), NULL),
  F13 = list(6, bf(y ~ a * x0 + c0, a ~ 1, c0 ~ s(x1), nl = TRUE),
             gaussian(), NULL),
  F14 = list(1, bf(y ~ s(x0) + s(x1) + s(x2) + s(x3)), gaussian(),
             y ~ s(x0) + s(x1) + s(x2) + s(x3)),
  F15 = list(2, bf(y ~ s(x, z, k = 30)), gaussian(), y ~ s(x, z, k = 30)),
  F16 = list(4, bf(y ~ fac + s(x2, by = fac) + s(x0)), gaussian(),
             y ~ fac + s(x2, by = fac) + s(x0)),
  F17 = list(1, bf(y ~ s(x0) + s(x2)), poisson(), y ~ s(x0) + s(x2)),
  F18 = list(6, bf(y ~ ti(x1) + ti(x2) + ti(x1, x2)), gaussian(),
             y ~ ti(x1) + ti(x2) + ti(x1, x2)))
rows <- list()
for (seed in seeds) {
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    d <- if (nm == "F17") gs(1, seed, re = FALSE, dist = "poisson", scale = 0.2) else
      if (cs[[1]] == 2) gs(2, seed, n = 300) else gs(cs[[1]], seed)
    nw <- 0L
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(withCallingHandlers(
      frm(cs[[2]], data = d, family = cs[[3]]),
      warning = function(x) {
        nw <<- nw + 1L
        invokeRestart("muffleWarning")
      }, message = function(m) invokeRestart("muffleMessage")),
      error = function(e) conditionMessage(e))
    el <- proc.time()[["elapsed"]] - t0
    if (is.character(fit)) {
      rows[[length(rows) + 1L]] <- list(seed = seed, case = nm, err = fit)
      next
    }
    fe <- suppressWarnings(fixef(fit))
    th <- fit$opt$par[names(fit$opt$par) == "theta"]
    ml <- NA_real_
    if (!is.null(cs[[4]])) {
      ml <- tryCatch({
        g <- mgcv::gam(cs[[4]], data = d, family = cs[[3]], method = "ML")
        -as.numeric(g$gcv.ubre)
      }, error = function(e) NA_real_)
    }
    rows[[length(rows) + 1L]] <- list(
      seed = seed, case = nm, err = NA_character_,
      conv = fit$opt$convergence, warns = nw, time = el,
      ll = as.numeric(logLik(fit)), ml = ml,
      se_ok = all(is.finite(fe[, "Est.Error"])),
      fe = fe[, "Estimate"], se = fe[, "Est.Error"], theta = th,
      fitted = as.numeric(fitted(fit)[, "Estimate"]))
  }
  cat("seed", seed, "done\n")
}
saveRDS(rows, a[4])
