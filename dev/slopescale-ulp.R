# Lane slopescale: do the new assertions of test-autoscale.R hold on
# every rounding? The data of each assertion moved by a few ulps,
# y * (1 + 2e-15 * k) (x for poisson, whose counts must stay integers),
# k = 1..20. Each row reports the quantity an assertion bounds, as a
# ratio to its bound (pass when below 1).
# Usage: Rscript slopescale-ulp.R [Rlib-arm-label]
lab <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", find.package("frmtmb"), "R", R.home(), "label", lab, "\n")
cat("BLAS", extSoftVersion()[["BLAS"]], "LAPACK", La_library(), "\n")
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
worst <- function(fs, f1, s) {
  slope <- function(nm) grepl("x$", nm)
  l1 <- as.numeric(logLik(f1))
  r <- c(ll = abs(as.numeric(logLik(fs)) - l1) / (1e-8 * abs(l1)))
  sds <- function(f, gn) {
    m <- VarCorr(f)[[gn]]$sd
    stats::setNames(m[, "Estimate"], rownames(m))
  }
  v <- 0
  for (gn in names(VarCorr(f1))) {
    a <- sds(f1, gn)
    b <- sds(fs, gn)
    b[slope(names(b))] <- b[slope(names(b))] * s
    v <- max(v, max(abs(b / a - 1)) / 1e-4)
  }
  e <- 0
  for (acc in list(ranef, coef)) {
    a1 <- acc(f1)
    as_ <- acc(fs)
    for (gn in names(a1)) {
      m1 <- as.matrix(a1[[gn]])
      ms <- as.matrix(as_[[gn]])
      ms <- sweep(ms, 2, ifelse(slope(colnames(ms)), s, 1), "*")
      e <- max(e, max(abs(ms - m1)) / max(abs(m1)) / 1e-4)
    }
  }
  y <- f1$frame$y[[1L]]
  ft <- max(abs(fitted(fs)[, "Estimate"] - fitted(f1)[, "Estimate"])) /
    stats::sd(y) / 1e-4
  c(r, sd = v, eff = e, fitted = ft)
}
rows <- list()
for (k in 1:20) {
  for (fam in c("gaussian", "poisson")) {
    for (form in c("0", "1")) {
      fo <- if (form == "0") y ~ 1 + (0 + x | g) else y ~ 1 + (1 + x | g)
      s <- 1e-4
      ds <- slope_data(2, s, fam)
      if (fam == "gaussian") ds$y <- ds$y * (1 + 2e-15 * k) else
        ds$x <- ds$x * (1 + 2e-15 * k)
      d1 <- ds
      d1$x <- ds$x / s
      fm <- if (fam == "gaussian") gaussian() else poisson()
      nw <- 0L
      wh <- function(w) {
        nw <<- nw + 1L
        invokeRestart("muffleWarning")
      }
      fs <- withCallingHandlers(frm(fo, family = fm, data = ds),
                                warning = wh)
      f1 <- withCallingHandlers(frm(fo, family = fm, data = d1),
                                warning = wh)
      rows[[length(rows) + 1L]] <- data.frame(
        k = k, case = paste(fam, form), t(worst(fs, f1, s)),
        autoscaled = isTRUE(fs$autoscaled), warnings = nw)
    }
  }
  # the guard absent: spread 1, default against autoscale = FALSE
  d <- slope_data(2, 1)
  d$y <- d$y * (1 + 2e-15 * k)
  fd <- frm(y ~ 1 + (1 + x | g), data = d)
  ff <- frm(y ~ 1 + (1 + x | g), data = d,
            control = frmtmb_control(autoscale = FALSE))
  # the verdict: a stop at the start is not called at zero
  d <- slope_data(2, 1e-5)
  d$y <- d$y * (1 + 2e-15 * k)
  conds <- character(0)
  f <- withCallingHandlers(
    frm(y ~ 1 + (0 + x | g), data = d,
        control = frmtmb_control(autoscale = FALSE)),
    message = function(m) {
      conds <<- c(conds, paste(class(m)[1], conditionMessage(m)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      conds <<- c(conds, paste(class(w)[1], conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  rows[[length(rows) + 1L]] <- data.frame(
    k = k, case = "guard+verdict", ll = 0, sd = 0, eff = 0, fitted = 0,
    autoscaled = identical(fd$opt$par, ff$opt$par),
    warnings = as.integer(
      !any(grepl("frmtmb_boundary_fit|is at zero", conds)) &&
        any(grepl("stopped short of the maximum", conds, fixed = TRUE))))
}
r <- do.call(rbind, rows)
options(width = 160)
print(r, digits = 3)
m <- r[r$case != "guard+verdict", ]
g <- r[r$case == "guard+verdict", ]
cat("\nslope cases:", nrow(m), "; worst ratio to bound: ll",
    signif(max(m$ll), 3), "sd", signif(max(m$sd), 3), "eff",
    signif(max(m$eff), 3), "fitted", signif(max(m$fitted), 3),
    "; not autoscaled", sum(!m$autoscaled), "; with warnings",
    sum(m$warnings > 0),
    "\nguard identical:", sum(g$autoscaled), "of", nrow(g),
    "; verdict short not zero:", sum(g$warnings), "of", nrow(g), "\n")
