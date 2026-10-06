# Reviewer of lane ordmix, final check: (1) the cost of the degenerate
# check at large n against the fit; (2) what a user sees on a probit fit
# the check leaves silent ("probit saturation"); (3) B4's identified
# count against the Hessian rank. Seed 20261007.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
gen <- function(n, c = 1, noise = rlogis, cuts = c(-1.5, 0, 1.5)) {
  x <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * c * x + 1, -0.8 * c * x - 1) + noise(n)
  data.frame(x = x, z = rnorm(n), y = 1L + rowSums(outer(lat, cuts, ">")))
}
el <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
cat("== (1) cost at large n\n")
for (n in c(5000, 50000)) {
  set.seed(20261007)
  d <- gen(n)
  t_fit <- el(f <- suppressWarnings(frm(bf(y ~ x),
                family = mixture(cumulative(), cumulative()), data = d)))
  t_chk <- el(for (i in 1:3) frmtmb:::mixture_ord_degeneracy(f, "y"))
  t_chk <- t_chk / 3
  t_fc <- el(suppressWarnings(frmtmb:::mixture_ord_fit_check(f, "y")))
  cat(sprintf("n=%d cum+cum: fit %.2f s, mixture_ord_degeneracy %.2f s, whole fit check %.2f s\n",
              n, t_fit, t_chk, t_fc))
  t_fit <- el(f2 <- suppressWarnings(frm(bf(y ~ cs(x)),
                family = mixture(sratio(), acat()), data = d)))
  t_chk <- el(frmtmb:::mixture_ord_degeneracy(f2, "y"))
  cat(sprintf("n=%d sratio+acat cs(x): fit %.2f s, mixture_ord_degeneracy %.2f s\n",
              n, t_fit, t_chk))
}

cat("== (2) probit saturation: what the user sees\n")
for (s in 1:20) {
  set.seed(s)
  d <- gen(500, c = 10, noise = rnorm, cuts = c(-7.5, 0, 7.5))
  w <- character(0)
  f <- tryCatch(withCallingHandlers(
    frm(bf(y ~ x), family = mixture(cumulative("probit"),
                                    cumulative("probit")), data = d),
    warning = function(e) {
      w <<- c(w, paste("[fit]", conditionMessage(e)))
      invokeRestart("muffleWarning")
    }), error = function(e) conditionMessage(e))
  if (is.character(f)) {
    cat(sprintf("seed %d ERROR %s\n", s, substr(f, 1, 90)))
    next
  }
  se <- withCallingHandlers(fixef(f)[, "Est.Error"], warning = function(e) {
    w <<- c(w, paste("[fixef]", conditionMessage(e)))
    invokeRestart("muffleWarning")
  })
  if (all(is.finite(se))) next
  invisible(withCallingHandlers(capture.output(summary(f)),
                                warning = function(e) {
    w <<- c(w, paste("[summary]", conditionMessage(e)))
    invokeRestart("muffleWarning")
  }))
  g <- f$obj$gr(f$opt$par)
  dg <- frmtmb:::mixture_ord_degeneracy(f, "y")
  cat(sprintf("seed %d: NaN SEs %d, gradient finite %s, reach %s, sharpen %s\n",
              s, sum(!is.finite(se)), all(is.finite(g)),
              paste(signif(dg$reach, 3), collapse = ","),
              paste(signif(dg$sharpen, 3), collapse = ",")))
  for (m in unique(w)) cat("   ", substr(m, 1, 170), "\n")
}

cat("== (3) B4: identified count against the Hessian rank, no predictor\n")
rank_of <- function(f) {
  H <- optimHess(f$opt$par, f$obj$fn, f$obj$gr)
  ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  sum(ev > 1e-6 * max(ev))
}
set.seed(20261007)
n <- 1500
g <- factor(sample(c("a", "b"), n, TRUE))
y4 <- sample(1:4, n, TRUE, prob = c(0.3, 0.2, 0.25, 0.25))
y3 <- sample(1:3, n, TRUE, prob = c(0.4, 0.35, 0.25))
d <- data.frame(g = g, y = ifelse(g == "a", y4, y3))
d$yh <- ifelse(runif(n) < 0.2, 0L, d$y)
d$yh2 <- ifelse(runif(n) < ifelse(d$g == "a", 0.1, 0.35), 0L, d$y)
cases <- list(
  plain = function() frm(bf(y ~ 1), family = mixture(cumulative(),
                                                      cumulative()),
                         data = d),
  gr = function() frm(bf(y | thres(gr = g) ~ 1),
                      family = mixture(cumulative(), cumulative()), data = d),
  hurdle = function() frm(bf(yh ~ 1),
                          family = mixture(hurdle_cumulative(),
                                           hurdle_cumulative()), data = d),
  hurdle_gr = function() frm(bf(yh2 | thres(gr = g) ~ 1),
                             family = mixture(hurdle_cumulative(),
                                              hurdle_cumulative()),
                             data = d))
for (nm in names(cases)) {
  w <- character(0)
  f <- withCallingHandlers(cases[[nm]](), warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  })
  fam <- f$spec$responses[[1]]$family
  nth <- fam$thres$nthres
  n_id <- sum(nth) + length(nth) * isTRUE(fam$extra_cat)
  cat(sprintf("%-10s nthres %s: parameters %d, lane's identified count %d, Hessian rank %d, warned %s\n",
              nm, paste(nth, collapse = "+"), length(f$opt$par), n_id,
              rank_of(f), any(grepl("no distributional parameter", w))))
}
