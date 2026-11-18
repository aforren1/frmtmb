# Lane slopetrig: do the new assertions of test-autoscale.R hold on
# every rounding? The four cases of "a slope column far from 0 against
# its spread engages the default" and the three of "a correlation run
# out at a short stop is not called a boundary", each over
# y * (1 + 2e-15 k) (x for the poisson counts), k = 0..20. Prints per
# case the number of draws on which each assertion held, and every
# warning text seen; then the gr(by =) plan test over
# x * (1 + 2e-15 k).
# Usage: Rscript slopetrig-ulp.R base|lane
arm <- commandArgs(TRUE)[1]
stopifnot(arm %in% c("base", "lane"))
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")),
    "lib:", dirname(find.package("frmtmb")), "BLAS:",
    extSoftVersion()[["BLAS"]], "LAPACK:", La_library(), "\n")
run <- function(expr) {
  w <- character(0)
  m <- character(0)
  f <- withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, paste(class(x)[1], conditionMessage(x)))
      invokeRestart("muffleMessage")
    })
  list(fit = f, w = w, m = m)
}
mk <- function(seed, s, c0) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  z <- stats::rnorm(300)
  eta <- 0.5 + stats::rnorm(20, 0, 0.5)[g] +
    stats::rnorm(20, 0, 0.5)[g] * z + 0.3 * z
  data.frame(g = g, x = c0 + s * z, y = eta + stats::rnorm(300))
}
slope_data <- function(seed, s, family = "gaussian") {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x0 <- stats::rnorm(300)
  eta <- 0.5 + stats::rnorm(20, 0, 0.5)[g] +
    stats::rnorm(20, 0, 0.5)[g] * x0
  y <- if (family == "gaussian") eta + stats::rnorm(300) else
    stats::rpois(300, exp(eta))
  data.frame(g = g, x = x0 * s, y = y)
}
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
allowed <- c("Standard errors are not available for 1 of 5 parameters. theta_3",
             "VarCorr() entries move along a direction",
             "Optimizer did not report convergence: false convergence (8)")
wseen <- character(0)
for (cs in list(list(2000, 3, 2, y ~ 1 + (1 + x | g)),
                list(2000, 3, 2, y ~ x + (1 + x | g)),
                list(5, 0.1, 1, y ~ 1 + (1 + x | g)),
                list(5, 0.1, 3, y ~ x + (1 + x | g)))) {
  ok <- integer(10)
  names(ok) <- c("engaged", "loglik", "sd", "re", "re_int", "coef",
                 "fitted", "warn_ok", "fixef_x", "se")
  for (k in 0:20) {
    d <- mk(cs[[3]], cs[[2]], cs[[1]])
    d$y <- d$y * (1 + 2e-15 * k)
    m <- mean(d$x)
    sx <- stats::sd(d$x)
    dz <- d
    dz$x <- (d$x - m) / sx
    rd <- run(frm(cs[[4]], data = d))
    rz <- run(frm(cs[[4]], data = dz))
    fd <- rd$fit
    fz <- rz$fit
    lz <- as.numeric(logLik(fz))
    rdm <- as.matrix(ranef(fd)$g)
    rzm <- as.matrix(ranef(fz)$g)
    w <- c(rd$w, rz$w)
    vd <- suppressWarnings(VarCorr(fd))$g$sd["x", "Estimate"]
    vz <- suppressWarnings(VarCorr(fz))$g$sd["x", "Estimate"]
    fx <- "x" %in% rownames(suppressWarnings(fixef(fd)))
    wseen <- c(wseen, w)
    ok <- ok + c(
      isTRUE(fd$autoscaled),
      lz - as.numeric(logLik(fd)) < 1e-8 * abs(lz),
      rel(vd * sx, vz) < 1e-4,
      rel(rdm[, "x"] * sx, rzm[, "x"]) < 1e-4,
      rel(rdm[, "Intercept"], rzm[, "Intercept"] - rzm[, "x"] * m / sx) < 1e-4,
      rel(as.matrix(coef(fd)$g)[, "x"] * sx,
          as.matrix(coef(fz)$g)[, "x"]) < 1e-4,
      rel(fitted(fd)[, "Estimate"], fitted(fz)[, "Estimate"]) < 1e-4,
      all(vapply(w, function(x) any(startsWith(x, allowed)), NA)),
      !fx || rel(fixef(fd)["x", c("Estimate", "Est.Error")] * sx,
                 fixef(fz)["x", c("Estimate", "Est.Error")]) < 1e-4,
      # the test's standard-error checks; the year with no fixed x is
      # left out there (counted as held)
      !(fx || cs[[1]] != 2000) || (
        (fx || rel(fixef(fd)["Intercept", "Est.Error"],
                   fixef(fz)["Intercept", "Est.Error"]) < 1e-3) &&
          rel(fitted(fd)[, "Est.Error"], fitted(fz)[, "Est.Error"]) < 1e-3))
  }
  cat("engage c0", cs[[1]], "seed", cs[[3]], deparse(cs[[4]]), ":",
      paste(names(ok), ok, sep = "=", collapse = " "), "of 21\n")
}
for (cs in list(list(4, 1e-4, "gaussian"), list(5, 1e-4, "gaussian"),
                list(5, 1e-5, "poisson"))) {
  n_nob <- n_short <- n_lost <- 0L
  for (k in 0:20) {
    d <- slope_data(cs[[1]], cs[[2]], cs[[3]])
    d$y <- d$y * (1 + 2e-15 * k)
    if (cs[[3]] == "poisson") {
      d$y <- round(d$y)
      d$x <- d$x * (1 + 2e-15 * k)
    }
    r <- run(frm(y ~ 1 + (1 + x | g), data = d,
                 family = if (cs[[3]] == "gaussian") gaussian() else
                   poisson(),
                 control = frmtmb_control(autoscale = FALSE)))
    lost <- frmtmb:::sdr_of(r$fit)$se_lost
    n_nob <- n_nob + !any(grepl("frmtmb_boundary_fit", r$m, fixed = TRUE))
    n_short <- n_short + any(grepl("stopped short of the maximum", r$w,
                                   fixed = TRUE))
    n_lost <- n_lost + identical(unname(lost[names(lost) == "theta_3"]),
                                 "shortcor")
  }
  cat("shortcor", paste(cs, collapse = " "), ": no boundary msg", n_nob,
      "short warning", n_short, "lost shortcor", n_lost, "of 21",
      if (cs[[3]] == "poisson") "(x, not the counts, moved by ulps)",
      "\n")
}
# "a gr(by =) block reads its ratio over the rows it covers": the plan
# only, x moved by ulps
n_ok <- 0L
for (k in 0:20) {
  set.seed(1)
  g <- factor(rep(1:20, each = 15))
  z <- stats::rnorm(300)
  d <- data.frame(g = g, h = factor(ifelse(as.integer(g) <= 10, "a", "b")),
                  x = (2000 + 3 * z) * (1 + 2e-15 * k))
  d$y <- 0.5 + stats::rnorm(20, 0, 0.5)[g] +
    (0.3 + stats::rnorm(20, 0, 0.5)[g]) * z + stats::rnorm(300)
  fr <- frm(y ~ 1 + (1 + x | gr(g, by = h)), data = d, dry_run = "frame")
  zz <- frmtmb:::autoscale_plan(fr)$y.mu$z
  want <- vapply(levels(d$h), function(l) {
    on <- d$h == l
    sqrt(mean(d$x[on]^2)) / stats::sd(d$x[on])
  }, 0)
  got <- vapply(zz, function(e) e$ratio, 0)
  n_ok <- n_ok + (length(zz) == 2L &&
                    max(abs(sort(got) / sort(want) - 1)) < 1e-12 &&
                    min(got) > frmtmb:::autoscale_slope_ratio &&
                    !is.null(frmtmb:::autoscale_decide(fr, frmtmb_control(),
                                                       FALSE)))
}
cat("grby plan: all assertions held on", n_ok, "of 21\n")
cat("warning texts seen in the engage cases:\n")
print(table(substr(wseen, 1, 90)))
