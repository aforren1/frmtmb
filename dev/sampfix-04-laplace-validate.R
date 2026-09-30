# Lane sampfix, script 04: laplace draws against full draws.
#
#   Rscript dev/sampfix-04-laplace-validate.R
#
# (a) What the lane computes from laplace draws, the population-level
#     posterior_epred(re_formula = NA) and conditional_effects(), must
#     agree with the same quantity from full draws of the same model up
#     to Monte Carlo error.
# (b) The two ways the refused quantity COULD have been computed, at
#     each draw's conditional modes b-hat(theta) or with b drawn from its
#     Laplace conditional N(b-hat(theta), H(theta)^-1), measured against
#     the full draws' in-sample posterior_epred(). This is the evidence
#     behind refusing instead (dev/sampfix-findings.md).
#
# Seeds: data 2024, sampler 11 (both arms), Laplace-conditional draws 5.

LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
.libPaths(c(LANE, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))

set.seed(2024L)
G <- 12L
nper <- 10L
dd <- data.frame(g = factor(rep(seq_len(G), each = nper)))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + rnorm(G, 0, 0.6)[as.integer(dd$g)] +
  rnorm(nrow(dd), 0, 0.8)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)

samp <- function(lap) {
  suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 4, iter = 2000, refresh = 0, seed = 11,
               laplace = lap)))
}
full <- samp(FALSE)
lap <- samp(TRUE)
cat("draws: full ", nrow(full$draws), " x ", ncol(full$draws), ", laplace ",
    nrow(lap$draws), " x ", ncol(lap$draws), "\n", sep = "")

cmp <- function(lab, a, b) {
  ma <- colMeans(a)
  mb <- colMeans(b)
  sb <- apply(b, 2L, sd)
  sa <- apply(a, 2L, sd)
  cat(sprintf("%-40s max |mean diff|/sd %.4f   sd ratio [%.4f, %.4f]\n",
              lab, max(abs(ma - mb) / sb), min(sa / sb), max(sa / sb)))
}

cat("\n(a) what laplace draws compute, against the full draws\n")
nd <- data.frame(x = seq(-2, 2, length.out = 9), g = factor(1, levels = 1:G))
ea <- posterior_epred(lap, newdata = nd, re_formula = NA)
eb <- posterior_epred(full, newdata = nd, re_formula = NA)
cmp("posterior_epred(newdata, re_formula = NA)", ea, eb)
ia <- posterior_epred(lap, re_formula = NA)
ib <- posterior_epred(full, re_formula = NA)
cmp("posterior_epred(re_formula = NA), 120 rows", ia, ib)
ca <- conditional_effects(lap)[[1L]]
cb <- conditional_effects(full)[[1L]]
cat(sprintf("%-40s max |estimate diff|/se %.4f   se ratio [%.4f, %.4f]\n",
            "conditional_effects(), 100 grid points",
            max(abs(ca$estimate__ - cb$estimate__) / cb$se__),
            min(ca$se__ / cb$se__), max(ca$se__ / cb$se__)))
cat("(b) population-level spread the Monte Carlo error allows: the full ",
    "draws' bulk ESS on b_Intercept is ",
    round(posterior::ess_bulk(full$draws[, "b_Intercept"])), ", laplace ",
    round(posterior::ess_bulk(lap$draws[, "b_Intercept"])), "\n", sep = "")

cat("\n(b) the refused in-sample posterior_epred(), two constructions\n")
obj <- fit$obj
rnd <- obj$env$random
keep <- obj$env$last.par.best
m_int <- frmtmb.sample:::draws_internal_matrix(lap)
outer <- setdiff(colnames(m_int), "lp__")
th <- m_int[, outer, drop = FALSE]
stopifnot(ncol(th) == length(obj$par))
set.seed(5L)
full_rows <- function(bmat) {
  # template order, then brms's scale and names
  tpl <- fit$frame[["par_template"]]
  lab <- frmtmb::brms_par_labels(fit, include_random = TRUE)
  M <- matrix(NA_real_, nrow(th), length(lab), dimnames = list(NULL, lab))
  pos <- 0L
  oc <- 0L
  for (cp in names(tpl)) {
    len <- length(tpl[[cp]])
    j <- pos + seq_len(len)
    if (cp == "b") M[, j] <- bmat else {
      M[, j] <- th[, oc + seq_len(len)]
      oc <- oc + len
    }
    pos <- pos + len
  }
  M <- frmtmb.sample:::draws_to_natural(M, fit)
  structure(list(stanfit = NULL, draws = cbind(M, lp__ = 0), fit = fit),
            class = "frmtmb_draws")
}
bhat <- matrix(NA_real_, nrow(th), length(rnd))
bsmp <- bhat
for (s in seq_len(nrow(th))) {
  obj$fn(th[s, ])
  par <- obj$env$last.par
  bh <- par[rnd]
  H <- as.matrix(obj$env$spHess(par, random = TRUE))
  R <- chol(H)
  bhat[s, ] <- bh
  bsmp[s, ] <- bh + backsolve(R, stats::rnorm(length(rnd)))
}
obj$env$last.par.best <- keep
ec <- posterior_epred(full_rows(bhat))
es <- posterior_epred(full_rows(bsmp))
ef <- posterior_epred(full)
cmp("at the conditional modes b-hat(theta)", ec, ef)
cmp("b from its Laplace conditional", es, ef)
cat("\nDONE\n")
