# Lane gpselimit: approximate gp() fits where frmtmb_gp_basis fires (length
# scale under a tenth of the spacing) and the ridge to the white limit is
# NOT flat. Three constructions:
#  prior: sweep row 3 (white without a prior) with a lognormal prior on
#    lscale centred on frac * spacing: the maximum is the prior's;
#  noicpt: the basis coefficients at their prior sd at frac * spacing,
#    alternating signs, `y ~ 0 + gp(x, k)`, noise 0.01;
#  short: sweep row 3 with the optimizer's iterations capped, so that it
#    stops on the way to the limit.
# Per fit: the length scale over the spacing, the objective along the
# ridge (log sdgp up D/2, log lscale down 1 per unit t) at t = -2, 2, 12
# minus the fit's, the SEs and lost reasons, and the conditions.
# Usage: Rscript gpselimit-notflat2.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]
gs_libs(arm)
ns <- asNamespace("frmtmb")
ridge <- function(fit, bk) {
  om <- ns$outer_par_map(fit)
  thpos <- which(om$comp == "theta")[bk$theta_idx]
  p <- fit$opt$par
  r <- numeric(length(p))
  r[thpos[1]] <- bk$gp_D / 2
  r[thpos[2]] <- -1
  f0 <- fit$obj$fn(p)
  v <- vapply(c(-2, 2, 12), function(t) fit$obj$fn(p + t * r) - f0, 0)
  fit$obj$fn(p)
  v
}
show <- function(lab, q) {
  fit <- q$fit
  if (inherits(fit, "error")) {
    cat("NOTFLAT2", lab, "ERROR", conditionMessage(fit), "\n")
    return(invisible())
  }
  bk <- Filter(function(b) b$covstruct == "hsgp", fit$frame$re_blocks)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  rg <- ridge(fit, bk)
  ss <- gs_ses(fit)
  tn <- grepl("^theta", names(ss$se))
  cat(sprintf(paste0("NOTFLAT2 %s %s | ll %.6f conv %d | ",
                     "ls/spacing %.4g sd %.4g | ridge f(t)-f0 t=-2 %.3g ",
                     "t=2 %.3g t=12 %.3g | se(theta) %s | lost [%s] | ",
                     "warn [%s] | msg [%s]\n"),
              arm, lab, as.numeric(logLik(fit)),
              fit$opt$convergence, exp(th[2]) / bk$gp_spacing, exp(th[1]),
              rg[1], rg[2], rg[3],
              paste(format(ss$se[tn], digits = 4), collapse = ","),
              paste(paste0(names(ss$lost), "=", ss$lost), collapse = ","),
              paste(unique(q$wclass), collapse = ";"),
              paste(unique(q$mclass), collapse = ";")))
  invisible(fit)
}
s <- gs_sweep_row(3)
f0 <- show("row3-free", gs_fit(bf(y ~ gp(x, k = s$k)), data = s$d,
                                family = gaussian()))
bk0 <- Filter(function(b) b$covstruct == "hsgp", f0$frame$re_blocks)[[1]]
sp <- bk0$gp_spacing
for (frac in c(0.01, 0.03, 0.05, 0.09)) {
  pr <- set_prior(sprintf("lognormal(%.10g, 0.05)", log(frac * sp)),
                  class = "lscale")
  show(sprintf("prior frac %.2f", frac),
       gs_fit(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian(),
              prior = pr))
}
for (it in c(2, 4, 6, 8, 10, 15)) {
  show(sprintf("short iter.max %d", it),
       gs_fit(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian(),
              control = frmtmb_control(optCtrl = list(iter.max = it),
                                       restarts = 0)))
}
n <- 200
set.seed(4100)
x <- sort(runif(n, 0, 10))
for (K in c(10, 30)) for (frac in c(0.05, 0.09)) {
  bk <- Filter(function(b) b$covstruct == "hsgp",
               suppressWarnings(suppressMessages(frm(
                 bf(y ~ 0 + gp(x, k = K)),
                 data = data.frame(x = x, y = rnorm(n)),
                 family = gaussian())))$frame$re_blocks)[[1]]
  ft <- suppressWarnings(suppressMessages(frm(
    bf(y ~ 0 + gp(x, k = K)), data = data.frame(x = x, y = rnorm(n)),
    family = gaussian())))
  lp <- ft$frame$linpreds[[ns$gp_block_lp(bk, ft$frame)]]
  Z <- as.matrix(ns$lp_Z(ft, lp))[, bk$c_idx, drop = FALSE]
  sds <- as.vector(ns$hsgp_sds(c(0, log(frac * bk$gp_spacing)),
                               bk$aux_omega, bk$gp_iso, ns$gp_cov_of(bk)))
  set.seed(4200 + K)
  b <- sds * rep(c(1, -1), length.out = length(sds))
  y <- as.vector(Z %*% b) + 0.01 * rnorm(n)
  show(sprintf("noicpt K %d frac %.2f", K, frac),
       gs_fit(bf(y ~ 0 + gp(x, k = K)), data = data.frame(x = x, y = y),
              family = gaussian()))
}
