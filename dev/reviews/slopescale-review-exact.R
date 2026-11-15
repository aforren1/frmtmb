# Review of lane slopescale, attack items 2 and 6: is the slope's
# rescaling an exact reparameterization, checked at a FIXED parameter
# point rather than at two optima? Per structure (the lane's
# slopescale-structs.R cases, seed 1 by default):
#   ds   = data with x * 1e-4 (the "original" model the user fits)
#   dsc  = ds with x divided by the plan's factor c (the scaled frame)
#   p_o  = the lane's reported optimum on ds; p_sc = autoscale_map(p_o,
#          plan, "scaled")
# Both models are evaluated at those points with optCtrl iter.max = 0
# (frm() returns the fit at its start), and every accessor is compared:
# objective, logLik, VarCorr sd and cor and Est.Error, ranef, coef,
# fitted, vcov, predict() draws seed for seed. Then the gap between the
# lane's two optima (x * 1e-4 against x) is set against the curvature:
# the spread-1 optimum mapped onto ds is evaluated there, and the change
# of objective is compared with the parameter distance.
# Usage: Rscript slopescale-review-exact.R [seeds]
a <- commandArgs(TRUE)
seeds <- if (length(a)) eval(parse(text = a[1])) else 1
only <- if (length(a) > 1) a[2] else NULL
.libPaths(c("C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", find.package("frmtmb"), "\n")
S <- 1e-4
mk <- function(seed) {
  set.seed(seed)
  ng <- 20; np <- 12
  g <- factor(rep(seq_len(ng), each = np))
  h <- factor(sample(seq_len(ng), ng * np, TRUE))
  a <- factor(rep(1:5, each = ng * np / 5))
  tf <- factor(rep(1:4, length.out = ng * np))
  x0 <- rnorm(ng * np)
  u0 <- rnorm(ng, 0, 0.5)[g]
  u1 <- rnorm(ng, 0, 0.6)[g]
  ua <- rnorm(5, 0, 0.4)[a]
  ut <- rnorm(ng * 4, 0, 0.5)[(as.integer(g) - 1) * 4 + as.integer(tf)]
  eta <- 0.5 + u0 + u1 * x0 + 0.5 * ua * x0 + ut * x0
  sig <- exp(-0.3 + rnorm(ng, 0, 0.4)[g] * x0)
  data.frame(y = rnorm(length(eta), eta, sig), g = g, h = h, a = a, tf = tf,
             x = x0)
}
A <- local({
  set.seed(99)
  M <- matrix(rnorm(20 * 30), 20)
  K <- crossprod(t(M)) / 30
  dimnames(K) <- list(as.character(1:20), as.character(1:20))
  K
})
cases <- list(
  us = bf(y ~ 1 + (1 + x | g)), us0 = bf(y ~ 1 + (0 + x | g)),
  dbar = bf(y ~ 1 + (1 + x || g)), cs = bf(y ~ 1 + cs(1 + x | g)),
  ar1 = bf(y ~ 1 + (1 | g) + ar1(0 + tf:x | g)),
  homdiag = bf(y ~ 1 + (1 | g) + homdiag(0 + x | g)),
  grnocor = bf(y ~ 1 + (1 + x | gr(g, cor = FALSE))),
  grcov = bf(y ~ 1 + (1 + x | gr(g, cov = A))),
  student = bf(y ~ 1 + (1 + x | gr(g, dist = "student"))),
  mm = bf(y ~ 1 + (1 + x | mm(g, h))),
  nested = bf(y ~ 1 + (1 + x | a/g)),
  dpar = bf(y ~ 1 + (1 | g), sigma ~ 1 + (0 + x | g)),
  idlink = bf(y ~ 1 + (1 + x | p | g), sigma ~ 1 + (1 | p | g)),
  xtwin = bf(y ~ x + (1 + x | g)))
if (!is.null(only)) cases <- cases[only]
quiet <- function(expr) suppressMessages(suppressWarnings(expr))
at0 <- function(f, d, st) {
  quiet(frm(f, data = d, data2 = list(A = A), start = st,
            control = frmtmb_control(autoscale = FALSE,
                                     optCtrl = list(iter.max = 0))))
}
rel <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_)
  max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-300))
}
slope <- function(nm) grepl("x$", nm)
vc_tab <- function(f, col) {
  v <- VarCorr(f)
  unlist(lapply(names(v), function(gn) {
    m <- v[[gn]]$sd
    stats::setNames(m[, col], paste0(gn, ":", rownames(m)))
  }))
}
cor_tab <- function(f) {
  v <- VarCorr(f)
  unlist(lapply(v, function(e) {
    if (is.null(e$cor)) return(NULL)
    as.numeric(e$cor[, "Estimate", ])
  }))
}
re_tab <- function(l, fac) {
  unlist(lapply(names(l), function(gn) {
    m <- as.matrix(l[[gn]])
    if (length(dim(m)) == 2L) {
      sweep(m, 2, ifelse(slope(colnames(m)), fac, 1), "*")
    } else {
      # brms layout: levels x stats x coefficients
      a3 <- l[[gn]]
      for (k in seq_len(dim(a3)[3])) {
        if (slope(dimnames(a3)[[3]][k])) a3[, , k] <- a3[, , k] * fac
      }
      a3
    }
  }))
}
rows <- list()
for (seed in seeds) for (nm in names(cases)) {
  f <- cases[[nm]]
  d1 <- mk(seed)
  ds <- d1
  ds$x <- d1$x * S
  fs <- quiet(frm(f, data = ds, data2 = list(A = A)))
  pl <- frmtmb:::autoscale_plan(fs$frame)
  zsc <- unlist(lapply(pl, function(p) vapply(p$z, `[[`, 0, "scale")))
  if (!length(zsc)) {
    cat(nm, seed, "no z plan\n")
    next
  }
  c0 <- zsc[1]
  one_c <- max(abs(zsc / c0 - 1))
  dsc <- ds
  dsc$x <- ds$x / c0
  po <- fs$estimates
  psc <- frmtmb:::autoscale_map(po, pl, "scaled")
  # the X plan (xtwin) centers and scales the fixed column; the data copy
  # then carries the z scaling only, so map only the z part
  if (nm == "xtwin") {
    plz <- lapply(pl, function(p) {
      p$cols <- integer(0)
      p
    })
    psc <- frmtmb:::autoscale_map(po, plz, "scaled")
    # dsc divides the fixed column by c0 too
    psc$beta[2] <- po$beta[2] * c0
  }
  fo0 <- at0(f, ds, po)
  fc0 <- at0(f, dsc, psc)
  r <- c(
    c_spread = one_c,
    obj = abs(fo0$opt$objective - fc0$opt$objective) /
      abs(fo0$opt$objective),
    vc_sd = rel(vc_tab(fc0, "Estimate") /
                  ifelse(slope(names(vc_tab(fc0, "Estimate"))), c0, 1),
                vc_tab(fo0, "Estimate")),
    vc_cor = {
      a1 <- cor_tab(fc0); b1 <- cor_tab(fo0)
      if (length(a1)) max(abs(a1 - b1)) else 0
    },
    vc_err = rel(vc_tab(fc0, "Est.Error") /
                   ifelse(slope(names(vc_tab(fc0, "Est.Error"))), c0, 1),
                 vc_tab(fo0, "Est.Error")),
    ranef = rel(re_tab(ranef(fc0), 1 / c0), re_tab(ranef(fo0), 1)),
    coef = rel(re_tab(coef(fc0), 1 / c0), re_tab(coef(fo0), 1)),
    fitted = rel(fitted(fc0)[, "Estimate"], fitted(fo0)[, "Estimate"]),
    vcov = rel(vcov(fc0), vcov(fo0)),
    pred = {
      set.seed(5); p1 <- predict(fo0)[, "Estimate"]
      set.seed(5); p2 <- predict(fc0)[, "Estimate"]
      max(abs(p1 - p2)) / sd(ds$y)
    })
  # the lane's comparison: the spread-1 optimum, mapped onto ds
  f1 <- quiet(frm(f, data = d1, data2 = list(A = A)))
  p1 <- f1$estimates
  # Z on ds is S times Z on d1: theta of a slope sd moves by log S, b / S
  pl1 <- lapply(pl, function(p) {
    p$cols <- integer(0)
    p$z <- lapply(p$z, function(z) {
      z$scale <- 1 / S
      z
    })
    p
  })
  # the slope's fixed effect (xtwin) is in units of 1/S too
  if (nm == "xtwin") p1$beta[2] <- p1$beta[2] / S
  p1o <- frmtmb:::autoscale_map(p1, pl1, "scaled")
  f1o <- at0(f, ds, p1o)
  th_s <- fs$opt$par
  th_1 <- f1o$opt$par
  r <- c(r, dobj_optima = f1o$opt$objective - fs$opt$objective,
         dpar_optima = max(abs(th_s - th_1)),
         hess_quad = {
           H <- tryCatch(stats::optimHess(fs$opt$par, fs$obj$fn, fs$obj$gr),
                         error = function(e) NULL)
           dv <- th_1 - th_s
           if (is.null(H)) NA else 0.5 * sum(dv * (H %*% dv))
         })
  rows[[length(rows) + 1L]] <- data.frame(case = nm, seed = seed,
                                          t(signif(r, 3)))
}
out <- do.call(rbind, rows)
options(width = 250)
print(out, row.names = FALSE)
ex <- c("obj", "vc_sd", "vc_cor", "vc_err", "ranef", "coef", "fitted", "vcov",
        "pred")
cat("\nat a fixed point, largest over cases:\n")
print(signif(apply(out[, ex], 2, max, na.rm = TRUE), 3))
