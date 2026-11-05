# Item 4: me() on new data at fixed parameter values, every shape, by reader.
# Data seed 7201, RNG seed 7202. Same values in both packages: the frm() estimate
# as every draw of a Fixed_param brmsfit and of frmtmb_draws.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7201)
n <- 300
d <- data.frame(z = rnorm(n), g = factor(rep(sprintf("g%02d", 1:30), each = 10)))
d$lat <- rnorm(n)
d$sx <- runif(n, 0.2, 0.6)
d$xo <- d$lat + rnorm(n, 0, d$sx)
d$lat2 <- 0.4 * d$lat + rnorm(n)
d$s2 <- 0.4
d$x2 <- d$lat2 + rnorm(n, 0, d$s2)
gl <- rnorm(30)[d$g]
d$sg <- 0.4
d$xg <- (gl + rnorm(30, 0, 0.4)[d$g])
d$y <- 1 + 1.2 * d$lat - 0.5 * d$z + rnorm(n, 0, 0.5)
d$y2 <- -0.5 + 0.9 * d$lat + rnorm(n, 0, 0.5)
d$yg <- 1 + 1.2 * gl - 0.5 * d$z + rnorm(n, 0, 0.5)
d$yc <- 1 + 1.2 * d$lat + 0.8 * d$lat2 + rnorm(n, 0, 0.5)
d$yi <- 1 + 1.2 * d$lat + 0.7 * d$lat * d$z + rnorm(n, 0, 0.5)
d$yp <- 1 + 0.9 * d$lat * d$lat2 + rnorm(n, 0, 0.5)
d$cnt <- rpois(n, exp(0.5 + 0.8 * d$lat))
d$bin <- rbinom(n, 1, plogis(0.3 + 1.2 * d$lat))
d$ys <- 0.5 * d$z + rnorm(n, 0, exp(-0.3 + 0.6 * d$lat))
d$ord <- cut(0.9 * d$lat + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
d$ynl <- 1 + 1.2 * d$lat + 0.3 * d$z + rnorm(n, 0, 0.5)
nd <- data.frame(xo = c(-1, 0.5, 1.5, 0.5), sx = c(0.3, 0.3, 0.8, 1.2), z = c(0, 1, -1, 0),
                 x2 = c(0, 1, -1, 0.5), s2 = c(0.4, 0.8, 0.4, 0.4),
                 xg = c(0.2, 0.2, 1, -1), sg = c(0.4, 0.4, 0.8, 0.6),
                 g = factor(c("g01", "g01", "g02", "g03"), levels = levels(d$g)))
shapes <- list(
  plain = list(f = "bf(y ~ me(xo, sx) + z)", fam = "gaussian()"),
  gr = list(f = "bf(yg ~ me(xg, sg, gr = g) + z)", fam = "gaussian()"),
  correlated = list(f = "bf(yc ~ me(xo, sx) + me(x2, s2))", fam = "gaussian()"),
  inter = list(f = "bf(yi ~ me(xo, sx) + me(xo, sx):z)", fam = "gaussian()"),
  me_me = list(f = "bf(yp ~ me(xo, sx):me(x2, s2))", fam = "gaussian()"),
  nl = list(f = "bf(ynl ~ a + b * z, a ~ me(xo, sx), b ~ 1, nl = TRUE)", fam = "gaussian()"),
  sigma = list(f = "bf(ys ~ z, sigma ~ me(xo, sx))", fam = "gaussian()"),
  log = list(f = "bf(cnt ~ me(xo, sx))", fam = "poisson()"),
  logit = list(f = "bf(bin ~ me(xo, sx))", fam = "bernoulli()"),
  ordinal = list(f = "bf(ord ~ me(xo, sx))", fam = "cumulative()"))
M <- 4000L
fmt <- function(v) paste(sprintf("%.4f", v), collapse = " ")
for (nm in names(shapes)) {
  sp <- shapes[[nm]]
  cat("\n=== shape", nm, sp$f, sp$fam, "\n")
  pr <- tryCatch(rev_pair(sp$f, sp$fam, d, ndraws = M), error = function(e) e)
  if (inherits(pr, "error")) { cat("PAIR ERROR:", conditionMessage(pr), "\n"); next }
  fit <- pr$fit; fx <- pr$fx
  ds <- rev_draws_at(fit, M)
  print(round(unlist(fit$estimates[c("beta", "betad")]), 4))
  set.seed(7202)
  if (nm != "ordinal") {
    bl <- quiet(posterior_linpred(fx, newdata = nd))
    be <- quiet(posterior_epred(fx, newdata = nd))
    fl <- try_num(frmtmb.sample::posterior_linpred(ds, newdata = nd))
    fe <- try_num(frmtmb.sample::posterior_epred(ds, newdata = nd))
    cat("linpred mean brms ", fmt(colMeans(bl)), "\n")
    if (is.matrix(fl)) cat("linpred mean draws", fmt(colMeans(fl)), "\n") else cat("draws ERR", attr(fl, "err"), "\n")
    cat("linpred sd   brms ", fmt(apply(bl, 2, sd)), "\n")
    if (is.matrix(fl)) cat("linpred sd   draws", fmt(apply(fl, 2, sd)), "\n")
    cat("epred mean brms ", fmt(colMeans(be)), "\n")
    if (is.matrix(fe)) cat("epred mean draws", fmt(colMeans(fe)), "\n") else cat("draws ERR", attr(fe, "err"), "\n")
    cat("epred sd   brms ", fmt(apply(be, 2, sd)), "\n")
    if (is.matrix(fe)) cat("epred sd   draws", fmt(apply(fe, 2, sd)), "\n")
    if (is.matrix(bl) && is.matrix(fl)) {
      cat("cor rows brms ", fmt(cor(bl)[1, ]), "\ncor rows draws", fmt(cor(fl)[1, ]), "\n")
    }
    pl <- try_num(frm_linpred(fit, newdata = nd, type = "link"))
    ps <- try_num(frm_linpred(fit, newdata = nd, type = "link", se.fit = TRUE))
    pb <- try_num(frm_lp_basis(fit, newdata = nd))
    fv <- try_num(fitted(fit, newdata = nd))
    cat("frm_linpred link     ", if (is.null(attr(pl, "err"))) fmt(pl) else attr(pl, "err"), "\n")
    cat("frm_linpred se       ", if (is.list(ps)) fmt(ps$se.fit) else attr(ps, "err"), "\n")
    cat("lp_basis extra sd    ", if (is.list(pb)) fmt(sqrt(pb$extra_var)) else attr(pb, "err"), "\n")
    cat("fitted Estimate      ", if (is.matrix(fv)) fmt(fv[, "Estimate"]) else attr(fv, "err"), "\n")
    cat("fitted Est.Error     ", if (is.matrix(fv)) fmt(fv[, "Est.Error"]) else "", "\n")
    cat("fitted Q2.5/Q97.5    ", if (is.matrix(fv)) paste(fmt(fv[, "Q2.5"]), "/", fmt(fv[, "Q97.5"])) else "", "\n")
    cat("brms epred Q2.5/Q97.5", fmt(apply(be, 2, quantile, .025)), "/", fmt(apply(be, 2, quantile, .975)), "\n")
  } else {
    be <- quiet(posterior_epred(fx, newdata = nd))
    fe <- try_num(frmtmb.sample::posterior_epred(ds, newdata = nd))
    fv <- try_num(fitted(fit, newdata = nd))
    cat("brms epred mean cat1", fmt(colMeans(be[, , 1])), " sd ", fmt(apply(be[, , 1], 2, sd)), "\n")
    if (length(dim(fe)) == 3L) cat("draws epred mean cat1", fmt(colMeans(fe[, , 1])), " sd ", fmt(apply(fe[, , 1], 2, sd)), "\n")
    if (length(dim(fv)) == 3L) {
      cat("fitted cat1 Estimate", fmt(fv[, "Estimate", 1]), " Est.Error ", fmt(fv[, "Est.Error", 1]), "\n")
    } else cat("fitted ERR", attr(fv, "err"), "\n")
  }
  bp <- quiet(posterior_predict(fx, newdata = nd))
  fp <- try_num(predict(fit, newdata = nd, summary = FALSE, ndraws = M, propagate_error = FALSE))
  fp2 <- try_num(frmtmb.sample::posterior_predict(ds, newdata = nd))
  cat("predict mean brms ", fmt(colMeans(bp)), " sd ", fmt(apply(bp, 2, sd)), "\n")
  if (is.matrix(fp)) cat("predict mean frmtmb", fmt(colMeans(fp)), " sd ", fmt(apply(fp, 2, sd)), "\n") else cat("predict ERR", attr(fp, "err"), "\n")
  if (is.matrix(fp2)) cat("predict mean draws ", fmt(colMeans(fp2)), " sd ", fmt(apply(fp2, 2, sd)), "\n") else cat("draws pp ERR", attr(fp2, "err"), "\n")
  sm <- try_num({
    s <- simulate(fit, nsim = M, seed = 3, newdata = nd)
    t(vapply(s, function(v) as.numeric(v), numeric(nrow(s))))
  })
  if (is.matrix(sm)) cat("simulate mean frmtmb", fmt(colMeans(sm)), " sd ", fmt(apply(sm, 2, sd)), "\n") else cat("simulate ERR", attr(sm, "err"), "\n")
}
