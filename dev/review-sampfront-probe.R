# Review of lane sampfront: attack the two-fill laplace probe/watch.
#   Rscript dev/review-sampfront-probe.R [base]
# Cases where the integrated values ARE read but the result is the same at
# fills 0 and 0.5772 (a linear combination of the fill that cancels, a flat
# function). Each case compares laplace draws with the full draws they come
# from; "SILENT WRONG" is a laplace answer that differs from the full one.
args <- commandArgs(TRUE)
base <- length(args) && identical(args[1], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", if (base) "base" else "lane", "frmtmb.sample",
    format(packageVersion("frmtmb.sample")), find.package("frmtmb.sample"),
    "\n")
lap_pair <- function(fit, n = 6L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  inner <- setdiff(lab, frmtmb::brms_par_labels(fit, include_random = FALSE))
  full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                    class = "frmtmb_draws")
  lap <- full
  lap$draws <- M[, setdiff(colnames(M), inner), drop = FALSE]
  list(full = full, lap = lap)
}
cmp <- function(lab, f) {
  cat("\n-- ", lab, "\n", sep = "")
  a <- tryCatch(suppressWarnings(f("lap")), error = function(e) e)
  b <- suppressWarnings(f("full"))
  if (inherits(a, "error")) {
    cat("  laplace: REFUSED:",
        substr(gsub("\\s+", " ", conditionMessage(a)), 1, 120), "\n")
    return(invisible())
  }
  d <- max(abs(unclass(a) - unclass(b)), na.rm = TRUE)
  cat("  laplace returned", length(a), "cells; max |laplace - full| =",
      format(d, digits = 6), if (d > 0) " SILENT WRONG" else " (equal)", "\n")
  cat("  laplace[1:3]:", format(head(as.numeric(a), 3), digits = 6),
      "\n  full[1:3]:   ", format(head(as.numeric(b), 3), digits = 6), "\n")
}
q <- function(...) suppressWarnings(suppressMessages(frm(...)))

set.seed(77)
G <- 8; n <- 12
dd <- data.frame(g = factor(rep(seq_len(G), each = n)))
dd$x <- rnorm(nrow(dd))
dd$f <- factor(rep(c("a", "b"), length.out = nrow(dd)))
u0 <- rnorm(G, 0, 0.8); u1 <- rnorm(G, 0, 0.8)
dd$y <- 0.5 + 0.4 * dd$x + u0[dd$g] + u1[dd$g] * dd$x + rnorm(nrow(dd), 0, 0.5)

# A. random intercept and slope, newdata at x = -1: u0 + u1 * x with a
# constant fill c is c * (1 - 1) = 0 at both fills, NA at fill NA
fa <- q(bf(y ~ x + (1 + x | g)), family = gaussian(), data = dd)
cat("par_template components:", names(fa$frame$par_template), "\n")
pa <- lap_pair(fa)
nda <- data.frame(x = -1, g = factor(c(1, 2), levels = 1:G))
cmp("A. epred (1 + x | g), newdata x = -1",
    function(w) posterior_epred(pa[[w]], newdata = nda))
cmp("A. linpred (1 + x | g), newdata x = -1",
    function(w) posterior_linpred(pa[[w]], newdata = nda))
cmp("A. predict (1 + x | g), newdata x = -1 (seeded)",
    function(w) { set.seed(5); posterior_predict(pa[[w]], newdata = nda) })
cmp("A control. epred (1 + x | g), newdata x = -0.5",
    function(w) posterior_epred(pa[[w]], newdata = transform(nda, x = -0.5)))

# B. sum-coded two-level factor with by-group slopes: the second level's
# design row is (1, -1)
op <- options(contrasts = c("contr.sum", "contr.poly"))
dd2 <- dd
dd2$y <- 0.5 + u0[dd2$g] + ifelse(dd2$f == "a", 1, -1) * u1[dd2$g] +
  rnorm(nrow(dd2), 0, 0.5)
fb <- q(bf(y ~ f + (1 + f | g)), family = gaussian(), data = dd2)
pb <- lap_pair(fb)
ndb <- data.frame(f = factor("b", levels = c("a", "b")),
                  g = factor(1:3, levels = 1:G))
cmp("B. epred (1 + f | g) sum-coded, newdata f = b",
    function(w) posterior_epred(pb[[w]], newdata = ndb))
options(op)

# C. uncorrelated (1 + x || g) at x = -1
fc <- q(bf(y ~ x + (1 + x || g)), family = gaussian(), data = dd)
pc <- lap_pair(fc)
cmp("C. epred (1 + x || g), newdata x = -1",
    function(w) posterior_epred(pc[[w]], newdata = nda))

# D. a breakpoint (hockey stick) nonlinear model with a group-level
# breakpoint: below it the term is (x - omega) * step(x - omega)
dd$yb <- 1 + 0.8 * pmax(dd$x - (0.2 + u0[dd$g] * 0.5), 0) +
  rnorm(nrow(dd), 0, 0.2)
fd <- tryCatch(q(bf(yb ~ a + bb * (x - om) * step(x - om), a ~ 1, bb ~ 1,
                    om ~ 1 + (1 | g), nl = TRUE), family = gaussian(),
                 data = dd,
                 prior = c(set_prior("normal(0, 2)", nlpar = "a"),
                           set_prior("normal(0, 2)", nlpar = "bb"),
                           set_prior("normal(0, 2)", nlpar = "om"))),
               error = function(e) e)
if (inherits(fd, "error")) {
  cat("\n-- D. step() nl fit failed:", conditionMessage(fd), "\n")
} else {
  pd <- lap_pair(fd)
  ndd <- data.frame(x = -3, g = factor(1:2, levels = 1:G))
  cmp("D. epred step() breakpoint, newdata x = -3",
      function(w) posterior_epred(pd[[w]], newdata = ndd))
  # a point just below the population breakpoint, in the group whose
  # breakpoint is lowest: the fills 0 and 0.5772 put it below the break
  # (term 0), the group's own value puts it above
  b_om <- unlist(fd$estimates$beta)
  cat("  beta:", paste(names(b_om), format(b_om, digits = 4)), "\n")
  om0 <- b_om[grep("om", names(b_om))][1]
  re <- fd$estimates$b
  cat("  b:", format(re, digits = 3), "\n")
  gmin <- which.min(re)
  ndd2 <- data.frame(x = om0 - 0.25, g = factor(gmin, levels = 1:G))
  cmp(sprintf("D2. epred step() breakpoint, x = om - 0.25 = %.3f, group %d",
              om0 - 0.25, gmin),
      function(w) posterior_epred(pd[[w]], newdata = ndd2))
}

# E. conditional_effects(re_formula = NULL) at x = -1 only
cmp("E. conditional_effects (1 + x | g), re_formula = NULL, x = -1",
    function(w) {
      ce <- conditional_effects(pa[[w]], effects = "x", re_formula = NULL,
                                int_conditions = list(x = -1),
                                conditions = data.frame(g = factor(1, levels = 1:G)))
      ce[[1]]$estimate__
    })
cat("\nDONE\n")
