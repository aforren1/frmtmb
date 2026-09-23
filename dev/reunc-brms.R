# Lane wt-reunc: a partial re_formula against brms on the same data,
# and on frmtmb.sample's draws.
#
# Three questions, each asked on the SAME data set:
#  1. What does brms keep? Its posterior_epred() under each re_formula
#     is checked against the formula it should equal, built by hand from
#     its own draws (b_Intercept + b_x x + the kept r_ columns).
#  2. Does frmtmb keep the same terms? The dropped contribution
#     fitted(NULL) - fitted(re_formula) is compared between frmtmb (at
#     the ML modes) and brms (the posterior mean). The two estimators
#     differ, so this is a measurement of closeness, reported beside the
#     size of the contribution itself; the base build keeps every term
#     and so reports a dropped contribution of exactly zero.
#  3. On draws (frmtmb.sample), posterior_epred() under each formula is
#     checked, draw by draw, against the same hand construction from the
#     draw's own b, and predict()/fitted() are the summaries of those.
#
# The brms fit is cached at dev/reunc-log/brms-partial-fit.rds so the
# base and the lane run read the SAME posterior. That file is 4 MB and
# is deleted after the two runs rather than committed; this script
# refits it in about two minutes.
#
#   Rscript dev/reunc-brms.R <lib>
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
cat("brms", format(packageVersion("brms")), "\n")
set.seed(3)
G <- 15; m <- 8
d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                h = factor(rep(1:5, times = G * m / 5)), x = rnorm(G * m))
ug <- rnorm(G, 0, 0.8); sg <- rnorm(G, 0, 0.4); uh <- rnorm(5, 0, 0.6)
d$y <- 1 + 0.5 * d$x + ug[d$g] + sg[d$g] * d$x + uh[d$h] + rnorm(G * m)

bfile <- file.path("dev/reunc-log", "brms-partial-fit.rds")
bfit <- if (file.exists(bfile)) readRDS(bfile) else {
  f <- brms::brm(y ~ x + (1 + x | g) + (1 | h), data = d, chains = 2,
                 iter = 3000, seed = 11, refresh = 0)
  saveRDS(f, bfile)
  f
}
fit <- frm(bf(y ~ x + (1 + x | g) + (1 | h)), data = d)
forms <- list("~(1 | g)" = ~(1 | g), "~(0 + x | g)" = ~(0 + x | g),
              "~(1 + x | g)" = ~(1 + x | g), "~(1 | h)" = ~(1 | h),
              "~(1 | g) + (1 | h)" = ~(1 | g) + (1 | h))
dr <- as.matrix(bfit)
rg <- function(lev, cf) dr[, sprintf("r_g[%s,%s]", lev, cf)]
rh <- function(lev) dr[, sprintf("r_h[%s,Intercept]", lev)]
hand <- function(ig, sl, ih) {
  out <- matrix(NA_real_, nrow(dr), nrow(d))
  for (i in seq_len(nrow(d))) {
    out[, i] <- dr[, "b_Intercept"] + dr[, "b_x"] * d$x[i] +
      ig * rg(d$g[i], "Intercept") + sl * rg(d$g[i], "x") * d$x[i] +
      ih * rh(d$h[i])
  }
  out
}
code <- list("~(1 | g)" = c(1, 0, 0), "~(0 + x | g)" = c(0, 1, 0),
             "~(1 + x | g)" = c(1, 1, 0), "~(1 | h)" = c(0, 0, 1),
             "~(1 | g) + (1 | h)" = c(1, 0, 1))
ep_full <- colMeans(brms::posterior_epred(bfit))
f_full <- fitted(fit)[, "Estimate"]
cat("\n1 and 2. brms against its own hand construction, and frmtmb's",
    "dropped contribution against brms's\n")
for (nm in names(forms)) {
  ep <- brms::posterior_epred(bfit, re_formula = forms[[nm]])
  cd <- code[[nm]]
  hd <- hand(cd[1], cd[2], cd[3])
  fr <- tryCatch(fitted(fit, re_formula = forms[[nm]])[, "Estimate"],
                 error = function(e) rep(NA_real_, nrow(d)))
  drop_b <- ep_full - colMeans(ep)
  drop_f <- f_full - fr
  cat(sprintf(paste0("%-20s brms vs hand max|d| %.2e | dropped: brms sd ",
                     "%.4f, frmtmb sd %.4f, max|brms - frmtmb| %.4f, ",
                     "cor %.4f\n"),
              nm, max(abs(ep - hd)), sd(drop_b), sd(drop_f),
              max(abs(drop_b - drop_f)),
              suppressWarnings(cor(drop_b, drop_f))))
}
ep_ns <- brms::posterior_epred(bfit, re_formula = ~(1 | nosuch))
ep_na <- brms::posterior_epred(bfit, re_formula = NA)
cat("brms ~(1 | nosuch) identical to re_formula = NA:",
    identical(ep_ns, ep_na), "\n")
ns <- tryCatch(fitted(fit, re_formula = ~(1 | nosuch)),
               error = function(e) e)
cat("frmtmb ~(1 | nosuch):", if (inherits(ns, "error")) {
  paste(class(ns)[1], conditionMessage(ns))
} else if (identical(ns, fitted(fit))) "identical to NULL" else "answered",
"\n")
nd <- d[c(1, 9, 17), c("x", "g")]
bnd <- tryCatch(colMeans(brms::posterior_epred(bfit, newdata = nd,
                                               re_formula = ~(1 | g))),
                error = function(e) conditionMessage(e))
fnd <- tryCatch(fitted(fit, newdata = nd, re_formula = ~(1 | g))[, 1],
                error = function(e) conditionMessage(e))
cat("newdata without h, ~(1 | g): brms", format(bnd), "\n",
    "                             frmtmb", format(fnd), "\n")

cat("\n3. frmtmb.sample draws\n")
set.seed(21)
ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 600, warmup = 300,
                                  seed = 21, refresh = 0))
fx <- as.matrix(ds$draws)
S <- nrow(fx)
for (nm in names(forms)) {
  ep <- posterior_epred(ds, re_formula = forms[[nm]])
  cd <- code[[nm]]
  # draw by draw, from the draw's own fit
  mx <- 0
  for (k in unique(round(seq(1, S, length.out = 25)))) {
    sh <- frmtmb.sample:::draws_fit_at(ds, k)
    fe <- sh$estimates[["beta"]]
    B <- sh$estimates[["b"]]
    bk <- sh$frame$re_blocks
    bg <- t(matrix(B[bk[[1]]$b_idx], 2))
    bh <- B[bk[[2]]$b_idx]
    hv <- fe[1] + fe[2] * d$x + cd[1] * bg[as.integer(d$g), 1] +
      cd[2] * bg[as.integer(d$g), 2] * d$x + cd[3] * bh[as.integer(d$h)]
    mx <- max(mx, abs(ep[k, ] - hv))
  }
  set.seed(5)
  pp <- posterior_predict(ds, re_formula = forms[[nm]])
  fi <- fitted(ds, re_formula = forms[[nm]])
  set.seed(5)
  pr <- predict(ds, re_formula = forms[[nm]])
  cat(sprintf(paste0("%-20s epred vs hand max|d| %.2e (25 draws) | ",
                     "fitted Est = colMeans(epred) %s | predict Est vs ",
                     "colMeans(pp) max|d| %.2e\n"),
              nm, mx, isTRUE(all.equal(unname(fi[, 1]),
                                       unname(colMeans(ep)))),
              max(abs(pr[, 1] - colMeans(pp)))))
}
ns2 <- tryCatch(posterior_epred(ds, re_formula = ~(1 | nosuch)),
                error = function(e) e)
cat("draws ~(1 | nosuch):", if (inherits(ns2, "error")) {
  paste(class(ns2)[1], substr(conditionMessage(ns2), 1, 80))
} else "answered", "\n")
