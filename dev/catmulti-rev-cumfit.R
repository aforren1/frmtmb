# Reviewer: the cumulative tied-threshold fix on FITS, in one process
# (rtmb-pitfalls item 21): each model is fitted with the lane's
# ord_cumulative_logpmf(), then with 0.71.0's body swapped into the
# namespace, and the estimates, logLik, and the tape's value and
# gradient at random points are compared bit for bit. Run under the
# stock R and under OpenBLAS 0.3.26 + LAPACK.
# Log dev/catmulti-rev-out/cumfit-<arm>.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
cat("BLAS:", extSoftVersion()[["BLAS"]], "| La_library:", La_library(), "\n")
ns <- asNamespace("frmtmb")
new_fun <- get("ord_cumulative_logpmf", ns)
old_fun <- eval(parse("dev/catmulti-rev-old-cum.R")[[1]][[3]])
environment(old_fun) <- ns
stopifnot(!identical(deparse(body(old_fun)), deparse(body(new_fun))))
use <- function(f) {
  unlockBinding("ord_cumulative_logpmf", ns)
  assign("ord_cumulative_logpmf", f, envir = ns)
  lockBinding("ord_cumulative_logpmf", ns)
}
set.seed(20261014)
n <- 400
dd <- data.frame(x = rnorm(n), z = rnorm(n),
                 g = factor(sample(letters[1:12], n, TRUE)),
                 grp = factor(sample(c("A", "B"), n, TRUE)))
u <- rnorm(12, 0, 0.6)
lat <- 0.8 * dd$x + u[as.integer(dd$g)] + rlogis(n)
dd$y <- as.integer(cut(lat, c(-Inf, -1.5, -0.3, 0.6, 1.8, Inf)))
dd$yo <- factor(dd$y, levels = 1:5, ordered = TRUE)
# hurdle_cumulative: zeros with probability plogis(-1 + 0.5 z)
dd$yh <- ifelse(runif(n) < plogis(-1 + 0.5 * dd$z), 0L, dd$y)
# grouped thresholds: group B's category 5 merged into 4
dd$yg <- ifelse(dd$grp == "B" & dd$y == 5, 4L, dd$y)
models <- list(
  logit_re = list(f = bf(y ~ x + (1 | g)), fam = cumulative()),
  probit = list(f = bf(y ~ x), fam = cumulative("probit")),
  cloglog_disc = list(f = bf(y ~ x, disc ~ z), fam = cumulative("cloglog")),
  thres_gr = list(f = bf(yg | thres(gr = grp) ~ x), fam = cumulative()),
  hurdle = list(f = bf(yh ~ x, hu ~ z), fam = hurdle_cumulative()),
  equidistant = list(f = bf(y ~ x), fam = cumulative(threshold = "equidistant"))
)
fitpair <- function(m) {
  out <- list()
  for (arm in c("new", "old")) {
    use(if (arm == "new") new_fun else old_fun)
    fit <- suppressMessages(suppressWarnings(frm(m$f, family = m$fam, data = dd)))
    set.seed(99)
    p0 <- fit$obj$env$last.par.best
    fixed <- fit$obj$par
    pts <- lapply(1:5, function(i) fixed + rnorm(length(fixed), 0, 0.3))
    out[[arm]] <- list(
      est = unlist(fit$estimates), ll = as.numeric(logLik(fit)),
      fn = vapply(pts, function(p) fit$obj$fn(p), 0),
      gr = lapply(pts, function(p) fit$obj$gr(p)),
      se = tryCatch(fixef(fit)[, "Est.Error"], error = function(e) NA))
  }
  use(new_fun)
  out
}
for (nm in names(models)) {
  r <- tryCatch(fitpair(models[[nm]]), error = function(e) {
    use(new_fun); cat(nm, "ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) next
  cat(sprintf(paste("%-13s estimates %s  logLik %s (%.10f)  fn at 5 points %s",
                    " gr at 5 points %s  SEs %s\n"), nm,
              identical(r$new$est, r$old$est), identical(r$new$ll, r$old$ll),
              r$new$ll, identical(r$new$fn, r$old$fn),
              identical(r$new$gr, r$old$gr), identical(r$new$se, r$old$se)))
}
cat("done\n")
# the swap must reach the fit: a perturbed body must move logLik
pert <- function(...) new_fun(...) * (1 + 1e-9)
environment(pert) <- environment()
use(pert)
fp <- suppressMessages(suppressWarnings(frm(models$probit$f, family = models$probit$fam, data = dd)))
use(new_fun)
fq <- suppressMessages(suppressWarnings(frm(models$probit$f, family = models$probit$fam, data = dd)))
cat(sprintf("swap reaches the fit: perturbed logLik %.12f vs %.12f (differ: %s)\n",
            as.numeric(logLik(fp)), as.numeric(logLik(fq)),
            !identical(as.numeric(logLik(fp)), as.numeric(logLik(fq)))))
