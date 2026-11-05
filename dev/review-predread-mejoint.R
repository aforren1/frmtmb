# Item 5: one me() term read by two responses or two dpars. Data seed 7701.
# What brms does on new data, what it does in sample with saved latents, what
# Stan's model says, and what the joint looks like either way.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7701)
n <- 400
d <- data.frame(lat = rnorm(n))
d$sx <- 0.5
d$xo <- d$lat + rnorm(n, 0, d$sx)
d$y1 <- 1 + 1.5 * d$lat + rnorm(n, 0, 0.4)
d$y2 <- -1 + 1.2 * d$lat + rnorm(n, 0, 0.4)
d$ys <- 0.5 + 1.0 * d$lat + rnorm(n, 0, exp(-0.8 + 0.5 * d$lat))
ftxt <- "mvbf(bf(y1 ~ me(xo, sx)), bf(y2 ~ me(xo, sx)), rescor = FALSE)"
fb <- rev_ev(ftxt, "brms")
cat("--- Stan's model: one latent for both responses?\n")
code <- brms::stancode(fb, data = d)
cat(grep("Xme|vector\\[Mme", strsplit(code, "\n")[[1L]], value = TRUE), sep = "\n")
pr <- rev_pair(ftxt, "gaussian()", d, ndraws = 4000)
fit <- pr$fit
fx <- pr$fx
ds <- rev_draws_at(fit, 4000L)
print(round(unlist(fit$estimates[c("beta", "betad")]), 4))
nd <- data.frame(xo = c(-1, 0, 1), sx = c(0.5, 1, 2))
set.seed(7702)
be <- posterior_epred(fx, newdata = nd)
bp <- posterior_predict(fx, newdata = nd)
cat("dim brms epred", dim(be), " dimnames3", dimnames(be)[[3]], "\n")
fmt <- function(v) paste(sprintf("%.4f", v), collapse = " ")
rowcor <- function(a) vapply(seq_len(dim(a)[2L]), function(j) cor(a[, j, 1], a[, j, 2]), 0)
cat("brms epred cor(y1,y2) by row:", fmt(rowcor(be)), "\n")
cat("brms predict cor(y1,y2) by row:", fmt(rowcor(bp)), "\n")
fe <- frmtmb.sample::posterior_epred(ds, newdata = nd)
fp <- frmtmb.sample::posterior_predict(ds, newdata = nd)
cat("draws epred cor by row:", fmt(rowcor(fe)), "\n")
cat("draws predict cor by row:", fmt(rowcor(fp)), "\n")
set.seed(2); pm <- predict(fit, newdata = nd, summary = FALSE, ndraws = 4000, propagate_error = FALSE)
cat("predict(fit) class/dim:", class(pm), dim(pm), "\n")
if (length(dim(pm)) == 3L) cat("predict(fit) cor by row:", fmt(rowcor(pm)), "\n")
b1 <- fit$estimates$beta[grep("mexosx", names(fit$estimates$beta))[1L]]
b2 <- fit$estimates$beta[grep("mexosx", names(fit$estimates$beta))[2L]]
s1 <- exp(fit$estimates$betad[1L])
s2 <- exp(fit$estimates$betad[2L])
cat("b1", b1, "b2", b2, "sigma1", s1, "sigma2", s2, "\n")
theo <- b1 * b2 * nd$sx^2 / sqrt((b1^2 * nd$sx^2 + s1^2) * (b2^2 * nd$sx^2 + s2^2))
cat("shared latent, predictive cor by row (theory):", fmt(theo), "\n")
cat("shared latent, epred cor by row (theory): 1\n")
cat("var(y1 + y2) brms       :", fmt(apply(bp[, , 1] + bp[, , 2], 2, var)), "\n")
cat("var(y1 + y2) shared theo:", fmt((b1 + b2)^2 * nd$sx^2 + s1^2 + s2^2), "\n")
cat("var(y1 - y2) brms       :", fmt(apply(bp[, , 1] - bp[, , 2], 2, var)), "\n")
cat("var(y1 - y2) shared theo:", fmt((b1 - b2)^2 * nd$sx^2 + s1^2 + s2^2), "\n")
# An interval for the difference of the two responses at one row, against the
# truth of the model: the model has one latent per row
cat("95% width of y1 - y2 at sdx = 2: indep", round(2 * 1.96 * sqrt((b1^2 + b2^2) * 4 + s1^2 + s2^2), 3),
    " shared", round(2 * 1.96 * sqrt((b1 - b2)^2 * 4 + s1^2 + s2^2), 3), "\n")
cat("\n--- the fit's own residual correlation in the training data\n")
cat("cor(y1, y2) in data, given xo, as a check of what the model implies:",
    round(cor(resid(lm(y1 ~ xo, d)), resid(lm(y2 ~ xo, d))), 4), "\n")
cat("posterior_predict at training rows as newdata, cor of y1 and y2 residuals over rows\n")
bt <- posterior_predict(fx, newdata = d[1:400, ])
ob1 <- d$y1 - colMeans(bt[, , 1])
ob2 <- d$y2 - colMeans(bt[, , 2])
cat("  observed residual correlation (y - E[y|x_obs]):", round(cor(ob1, ob2), 4), "\n")
cat("  brms replicate residual correlation (median over draws):",
    round(median(vapply(1:200, function(s) cor(bt[s, , 1] - colMeans(bt[, , 1]), bt[s, , 2] - colMeans(bt[, , 2])), 0)), 4), "\n")
th <- b1 * b2 * 0.25 / sqrt((b1^2 * 0.25 + s1^2) * (b2^2 * 0.25 + s2^2))
cat("  shared latent replicate correlation (theory at sdx 0.5):", round(th, 4), "\n")
cat("\n--- saved latents, in sample: one latent, so brms reads it for both\n")
bfs <- quiet(brms::brm(fb, data = d, chains = 1, iter = 500, seed = 3, refresh = 0, silent = 2,
                       save_pars = brms::save_pars(latent = TRUE)))
cat("variables with Xme:", paste(head(grep("^Xme", brms::variables(bfs), value = TRUE), 3), collapse = " "), "... count",
    length(grep("^Xme", brms::variables(bfs))), " n =", n, "\n")
a <- posterior_epred(bfs)
cat("in-sample, saved latents, cor(y1 mean, y2 mean) over draws, first 5 rows:",
    fmt(vapply(1:5, function(j) cor(a[, j, 1], a[, j, 2]), 0)), "\n")
b <- posterior_epred(bfs, newdata = d[1:5, ])
cat("same rows as newdata, same fit:",
    fmt(vapply(1:5, function(j) cor(b[, j, 1], b[, j, 2]), 0)), "\n")
cat("\n--- two dpars\n")
ptx <- "bf(ys ~ me(xo, sx), sigma ~ me(xo, sx))"
P2 <- rev_pair(ptx, "gaussian()", d, ndraws = 4000)
ds2 <- rev_draws_at(P2$fit, 4000L)
print(round(unlist(P2$fit$estimates[c("beta", "betad")]), 4))
nd2 <- data.frame(xo = c(-1, 0, 1), sx = c(0.5, 1, 2))
bmu <- posterior_epred(P2$fx, newdata = nd2)
bsig <- posterior_linpred(P2$fx, newdata = nd2, dpar = "sigma")
cat("brms cor(mu, log sigma), separately drawn calls:", fmt(vapply(1:3, function(j) cor(bmu[, j], bsig[, j]), 0)), "\n")
ppb <- posterior_predict(P2$fx, newdata = nd2)
ppf <- frmtmb.sample::posterior_predict(ds2, newdata = nd2)
cat("predict var brms :", fmt(apply(ppb, 2, var)), "\n")
cat("predict var draws:", fmt(apply(ppf, 2, var)), "\n")
kurt <- function(v) mean((v - mean(v))^4) / var(v)^2
cat("predict kurtosis brms :", fmt(apply(ppb, 2, kurt)), "\n")
cat("predict kurtosis draws:", fmt(apply(ppf, 2, kurt)), "\n")
bb <- P2$fit$estimates$beta
bd <- P2$fit$estimates$betad
mu_slope <- bb[grep("mexosx", names(bb))]
sg0 <- bd[1L]
sg_slope <- bd[grep("mexosx", names(bd))]
sim_pp <- function(x, s, shared, M = 200000L) {
  z1 <- rnorm(M)
  z2 <- if (shared) z1 else rnorm(M)
  m <- bb[1L] + mu_slope * (x + s * z1)
  sg <- exp(sg0 + sg_slope * (x + s * z2))
  m + rnorm(M, 0, sg)
}
set.seed(7703)
for (shared in c(FALSE, TRUE)) {
  v <- lapply(seq_len(3), function(i) sim_pp(nd2$xo[i], nd2$sx[i], shared))
  cat(if (shared) "shared" else "indep ", " predictive var", fmt(sapply(v, var)),
      " q99.5-q.5 width", fmt(sapply(v, function(q) diff(quantile(q, c(.005, .995))))), " skew",
      fmt(sapply(v, function(q) mean((q - mean(q))^3) / sd(q)^3)), "\n")
}
