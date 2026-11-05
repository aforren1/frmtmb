# What brms 2.23.0 does with me() on new data, RUN at fixed parameter values
# (Fixed_param through brms_fixed_fit(), 2000 equal draws, data seed 41):
# the spread over draws, whether two responses or two dpars share a latent,
# me(gr =) with new and missing levels, correlated me(), and the refusals.
#   Rscript dev/predread-me-brms-probe.R > dev/predread-log/me-brms-probe.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
  library(testthat)
})
E <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = E)
}
brms_fixed_fit <- E$brms_fixed_fit
tr <- function(lab, e) {
  cat("\n---", lab, "\n")
  r <- tryCatch(e, error = function(err) {
    cat("ERROR:", conditionMessage(err), "\n")
    invisible(NULL)
  })
  invisible(r)
}
set.seed(41)
n <- 120
g <- rep(1:12, each = 10)
d <- data.frame(g = factor(g), z = rnorm(n))
lat <- rnorm(12)[g]
d$xg <- lat + rnorm(12, 0, 0.3)[g]
d$sg <- 0.3
d$x <- rnorm(n)
d$sx <- 0.4
d$w <- rnorm(n)
d$sw <- 0.25
d$y <- 1 + 0.8 * d$x + 0.5 * d$w + 0.6 * lat + 0.3 * d$z + rnorm(n, 0, 0.5)
d$y2 <- -0.5 + 0.7 * d$x + rnorm(n, 0, 0.6)
nd <- data.frame(x = c(-1, 0, 1), sx = c(0.2, 0.4, 0.8), w = 0.5, sw = 0.25,
                 z = 0, xg = c(0.5, 0.5, -1), sg = c(0.3, 0.3, 0.3),
                 g = factor(c("1", "1", "2"), levels = levels(d$g)))

# one me()
fit1 <- frm(bf(y ~ me(x, sx) + z), data = d, family = gaussian())
f1 <- brms_fixed_fit(brms::bf(y ~ me(x, sx) + z), stats::gaussian(), d, fit1,
                     ndraws = 2000)
b <- brms::fixef(f1)[, "Estimate"]
print(b)
ep <- brms::posterior_epred(f1, newdata = nd)
cat("epred mean by row:", colMeans(ep), "\n")
cat("epred sd by row:", apply(ep, 2, sd), "\n")
cat("|b| sdx:", abs(b[["mexsx"]]) * nd$sx, "\n")
pl <- brms::posterior_linpred(f1, newdata = nd)
cat("linpred sd by row:", apply(pl, 2, sd), "\n")
cat("rows independent? cor(row1, row2):", cor(ep[, 1], ep[, 2]), "\n")
pp <- brms::posterior_predict(f1, newdata = nd)
cat("predict sd by row:", apply(pp, 2, sd), "\n")
fi <- fitted(f1, newdata = nd)
print(fi)
tr("sdx = 0 in new data", {
  nd0 <- nd
  nd0$sx <- 0
  print(apply(brms::posterior_epred(f1, newdata = nd0), 2, sd))
})
tr("no sdx column in new data", {
  print(apply(brms::posterior_epred(f1, newdata = nd[c("x", "z")]), 2, sd))
})
tr("NA in x", {
  nd0 <- nd
  nd0$x[1] <- NA
  print(apply(brms::posterior_epred(f1, newdata = nd0), 2, sd))
})
tr("in sample: spread of the fixed draws (latent draws are fixed)", {
  print(apply(brms::posterior_epred(f1)[, 1:3], 2, sd))
})
tr("emmeans with sx not asked for", {
  print(as.data.frame(emmeans::emmeans(f1, ~ z, at = list(z = c(0, 1)))))
})

# two responses reading the same me(): does a row share its latent?
fit2 <- frm(mvbf(bf(y ~ me(x, sx)), bf(y2 ~ me(x, sx)), rescor = FALSE),
            data = d, family = gaussian())
f2 <- tr("mv fixed", brms_fixed_fit(
  brms::mvbf(brms::bf(y ~ me(x, sx)), brms::bf(y2 ~ me(x, sx)),
             rescor = FALSE), stats::gaussian(), d, fit2, ndraws = 2000))
tr("mv epred", {
  e2 <- brms::posterior_epred(f2, newdata = nd)
  print(dim(e2))
  cat("cor(y, y2) row 1..3 over draws:",
      sapply(1:3, function(i) cor(e2[, i, 1], e2[, i, 2])), "\n")
  cat("sd y:", apply(e2[, , 1], 2, sd), " sd y2:", apply(e2[, , 2], 2, sd),
      "\n")
})

# mu and sigma reading the same me()
fit3 <- frm(bf(y ~ me(x, sx), sigma ~ me(x, sx)), data = d, family = gaussian())
f3 <- tr("dpar fixed", brms_fixed_fit(
  brms::bf(y ~ me(x, sx), sigma ~ me(x, sx)), stats::gaussian(), d, fit3,
  ndraws = 2000))
tr("dpar: mu and sigma", {
  mu <- brms::posterior_linpred(f3, newdata = nd, dpar = "mu")
  sg <- brms::posterior_linpred(f3, newdata = nd, dpar = "sigma")
  cat("sd mu:", apply(mu, 2, sd), " sd log sigma:", apply(sg, 2, sd), "\n")
  nd$y <- 0
  pr <- brms::prepare_predictions(f3, newdata = nd)
  mu <- brms:::get_dpar(pr, "mu")
  sg <- brms:::get_dpar(pr, "sigma")
  cat("one prep, cor(mu, log sigma):",
      sapply(1:3, function(i) cor(mu[, i], log(sg[, i]))), "\n")
})

# me(gr = )
fit4 <- frm(bf(y ~ me(xg, sg, gr = g) + z), data = d, family = gaussian())
f4 <- tr("gr fixed", brms_fixed_fit(
  brms::bf(y ~ me(xg, sg, gr = g) + z), stats::gaussian(), d, fit4,
  ndraws = 2000))
tr("gr: rows of one level share a latent", {
  e4 <- brms::posterior_epred(f4, newdata = nd)
  cat("sd by row:", apply(e4, 2, sd), "\n")
  cat("cor(row 1, row 2) same level:", cor(e4[, 1], e4[, 2]), "\n")
  cat("cor(row 1, row 3) other level:", cor(e4[, 1], e4[, 3]), "\n")
})
tr("gr: new level", {
  nd2 <- nd
  nd2$g <- as.character(nd2$g)
  nd2$g[3] <- "new"
  e4 <- brms::posterior_epred(f4, newdata = nd2, allow_new_levels = TRUE)
  cat("sd by row:", apply(e4, 2, sd), "\n")
})
tr("gr: new level, not allowed", {
  nd2 <- nd
  nd2$g <- as.character(nd2$g)
  nd2$g[3] <- "new"
  e4 <- brms::posterior_epred(f4, newdata = nd2)
  cat("sd by row:", apply(e4, 2, sd), "\n")
})
tr("gr: NA level", {
  nd2 <- nd
  nd2$g[c(1, 2)] <- NA
  e4 <- brms::posterior_epred(f4, newdata = nd2, allow_new_levels = TRUE)
  cat("sd by row:", apply(e4, 2, sd), " cor(1, 2):",
      cor(e4[, 1], e4[, 2]), "\n")
})
tr("gr: x varies in a level", {
  nd2 <- nd
  nd2$xg[2] <- 9
  print(apply(brms::posterior_epred(f4, newdata = nd2), 2, sd))
})
tr("gr: no g column", {
  print(apply(brms::posterior_epred(f4, newdata = nd[c("xg", "sg", "z")]), 2,
              sd))
})

# correlated me(): corme is not read on new data
fit5 <- frm(bf(y ~ me(x, sx) + me(w, sw)), data = d, family = gaussian())
f5 <- tr("cor fixed", brms_fixed_fit(
  brms::bf(y ~ me(x, sx) + me(w, sw)), stats::gaussian(), d, fit5,
  ndraws = 2000))
tr("cor: rows' latents independent per call", {
  nd$y <- 0
  pr <- brms::prepare_predictions(f5, newdata = nd)
  sp <- pr$dpars$mu$sp
  cat("names:", names(sp$Xme), "\n")
  cat("cor(Xme_x, Xme_w) by row:",
      sapply(1:3, function(i) cor(sp$Xme[[1]][, i], sp$Xme[[2]][, i])), "\n")
  cat("mean Xme_x:", colMeans(sp$Xme[[1]]), " sd:",
      apply(sp$Xme[[1]], 2, sd), "\n")
  cors <- grep("^cor", brms::variables(f5), value = TRUE)
  cat("corme in the draws:", cors,
      unlist(brms::as_draws_df(f5, variable = cors)[1, cors]), "\n")
})
# an interaction reads one latent in every column of a call
fit6 <- frm(bf(y ~ me(x, sx) + me(x, sx):z), data = d, family = gaussian())
f6 <- tr("interaction fixed", brms_fixed_fit(
  brms::bf(y ~ me(x, sx) + me(x, sx):z), stats::gaussian(), d, fit6,
  ndraws = 2000))
tr("interaction", {
  nd6 <- nd
  nd6$z <- c(1, 1, 1)
  e6 <- brms::posterior_epred(f6, newdata = nd6)
  bb <- brms::fixef(f6)[, "Estimate"]
  print(bb)
  cat("sd by row:", apply(e6, 2, sd), " one latent says |b1 + b2| sdx:",
      abs(bb[["mexsx"]] + bb[["mex:z"]]) * nd6$sx, "\n")
})
