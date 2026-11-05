# Item 4, ML analog: does the fit's reading of a new me() row match the truth?
#   Rscript dev/review-predread-sesim.R <seed_from> <seed_to> <out>
# Per dataset (seed): simulate from the me() model, frm() it, read a grid of new
# rows. Three targets:
#   T1 eta at the true latent x* where x* ~ N(x, sdx) is brms's draw (the claim)
#   T2 eta at x* where x_obs | x* ~ N(x*, sdx), x* ~ N(0, 1): the model's own
#      conditional law of a latent given what was observed
#   T3 E[y] at T1 on the response scale, for the log link
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
s1 <- as.integer(args[1])
s2 <- as.integer(args[2])
outf <- args[3]
nd <- data.frame(xo = c(-1.5, 0, 1.5, 0, 1.5), sx = c(0.3, 0.3, 0.3, 1, 1))
rows <- list()
for (s in s1:s2) {
  set.seed(s)
  n <- 200
  lat <- rnorm(n)
  sx <- 0.4
  xo <- lat + rnorm(n, 0, sx)
  d <- data.frame(xo = xo, sx = sx)
  d$y <- 1 + 1.0 * lat + rnorm(n, 0, 0.5)
  d$cnt <- rpois(n, exp(0.5 + 0.7 * lat))
  fit <- tryCatch(quiet(frm(bf(y ~ me(xo, sx)), family = gaussian(), data = d)), error = function(e) NULL)
  fitp <- tryCatch(quiet(frm(bf(cnt ~ me(xo, sx)), family = poisson(), data = d)), error = function(e) NULL)
  if (is.null(fit) || is.null(fitp)) next
  p <- frm_linpred(fit, newdata = nd, se.fit = TRUE)
  pp <- frm_linpred(fitp, newdata = nd, se.fit = TRUE)
  fv <- fitted(fitp, newdata = nd)
  # brms's draw of the latent
  z <- rnorm(nrow(nd))
  xd <- nd$xo + nd$sx * z
  # the model's conditional law of the latent, prior N(0, 1)
  v <- nd$sx^2 / (1 + nd$sx^2)
  xc <- nd$xo * 1 / (1 + nd$sx^2) + sqrt(v) * z
  T1 <- 1 + xd
  T2 <- 1 + xc
  Tp1 <- 0.5 + 0.7 * xd
  Tr <- exp(Tp1)
  # a predictive replicate at the draw, for predict()'s interval
  yrep <- T1 + rnorm(nrow(nd), 0, 0.5)
  set.seed(s + 1e5); pr <- predict(fit, newdata = nd, ndraws = 400)
  rows[[length(rows) + 1L]] <- data.frame(
    seed = s, row = seq_len(nrow(nd)), xo = nd$xo, sx = nd$sx,
    eta = p$fit, se = p$se.fit,
    c1 = abs(T1 - p$fit) <= 1.96 * p$se.fit, c2 = abs(T2 - p$fit) <= 1.96 * p$se.fit,
    err1 = T1 - p$fit, err2 = T2 - p$fit,
    etap = pp$fit, sep = pp$se.fit, cp1 = abs(Tp1 - pp$fit) <= 1.96 * pp$se.fit,
    f_est = fv[, "Estimate"], f_err = fv[, "Est.Error"], f_lo = fv[, "Q2.5"], f_hi = fv[, "Q97.5"],
    resp_cover = Tr >= fv[, "Q2.5"] & Tr <= fv[, "Q97.5"],
    resp_true_mean = exp(0.5 + 0.7 * nd$xo + 0.5 * 0.49 * nd$sx^2),
    pr_lo = pr[, "Q2.5"], pr_hi = pr[, "Q97.5"],
    pr_cover = yrep >= pr[, "Q2.5"] & yrep <= pr[, "Q97.5"],
    pr_sd = pr[, "Est.Error"],
    bhat = fit$estimates$beta[2L])
}
res <- do.call(rbind, rows)
saveRDS(res, outf)
cat("seeds", length(unique(res$seed)), "\n")
