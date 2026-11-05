# Short check of fix round 2, item 2: compat me() x marginaleffects and NEWS say "A slope,
# contrast or test of a coefficient is right". Measure it on an identity and a log link,
# per row and averaged, against draws of the coefficients with the new row's latent value
# drawn from N(xo, sx), which is what brms's draws hold. Seed 7431. 4000 draws.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7431)
n <- 300
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.5)
d$xo <- d$lat + rnorm(n, 0, 0.5)
d$y <- 1 + 0.8 * d$lat - 0.5 * d$z + rnorm(n, 0, 0.7)
d$cnt <- rpois(n, exp(0.5 + 0.6 * d$lat - 0.5 * d$z))
nd <- data.frame(xo = c(0, 1, 2), z = 0, sx = 0.5)
q <- function(e) suppressWarnings(suppressMessages(e))
draws <- function(fit, nd, inv, S = 4000L) {
  cf <- fixef(fit)[, "Estimate"]
  V <- as.matrix(stats::vcov(fit))[names(cf), names(cf)]
  B <- MASS::mvrnorm(S, cf, V)
  a <- B[, "Intercept"]; z <- B[, "z"]; b <- B[, grep("xo", colnames(B))]
  rows <- vapply(seq_len(nrow(nd)), function(i) {
    l <- rnorm(S, nd$xo[i], nd$sx[i])
    z * inv(a + b * l + z * nd$z[i])
  }, numeric(S))
  rows
}
for (nm in c("gaussian", "poisson")) {
  fit <- quiet(frm(bf(if (nm == "gaussian") y ~ me(xo, sx) + z else cnt ~ me(xo, sx) + z),
                   family = if (nm == "gaussian") gaussian() else poisson(), data = d))
  inv <- if (nm == "gaussian") function(e) rep(1, length(e)) else exp
  r <- draws(fit, nd, inv)
  s <- q(marginaleffects::slopes(fit, newdata = nd, variables = "z"))
  av <- q(marginaleffects::avg_slopes(fit, newdata = nd, variables = "z"))
  cat(sprintf("%-9s per row  marginaleffects est %s se %s\n", nm,
              paste(signif(s$estimate, 4), collapse = " "), paste(signif(s$std.error, 4), collapse = " ")))
  cat(sprintf("%-9s per row  draws           mean %s sd %s\n", nm,
              paste(signif(colMeans(r), 4), collapse = " "), paste(signif(apply(r, 2L, sd), 4), collapse = " ")))
  cat(sprintf("%-9s average  marginaleffects est %.4f se %.4f   draws mean %.4f sd %.4f\n", nm,
              av$estimate, av$std.error, mean(rowMeans(r)), sd(rowMeans(r))))
}
