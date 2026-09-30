# Reviewer no-regression run (lane ceplot): conditional_effects(),
# fitted() and predict() on 30 models that use no new option, saved per
# arm and compared with identical() by dev/ceplot-rev-noreg-compare.R.
#   Rscript dev/ceplot-rev-noreg.R lane|base
# Seeds: data 21 + case index; every stochastic call is preceded by
# set.seed(1) (or takes seed = 1).
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib",
          "C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
mk <- function(i, n = 200) {
  set.seed(21 + i)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  f = factor(sample(c("a", "b", "c"), n, TRUE)),
                  g = factor(rep(1:10, length.out = n)),
                  h = factor(rep(1:5, each = 2, length.out = n)),
                  o = sample(1:4, n, TRUE), sdx = runif(n, 0.1, 0.3),
                  tr = sample(5:10, n, TRUE))
  d$g1 <- factor(sample(1:10, n, TRUE), levels = 1:10)
  d$g2 <- factor(sample(1:10, n, TRUE), levels = 1:10)
  d$fg <- factor(ifelse(as.integer(d$g) <= 5, "a", "b"))
  u <- rnorm(10)
  eta <- 0.3 + 0.5 * d$x - 0.3 * d$z + u[d$g]
  d$y <- rnorm(n, eta, exp(0.1 * d$z))
  d$yp <- abs(d$y) + 0.2
  d$yc <- rpois(n, exp(0.2 + 0.3 * d$x + 0.3 * u[d$g]))
  d$yb <- rbinom(n, d$tr, plogis(0.3 * d$x + 0.3 * u[d$g]))
  d$yo <- factor(cut(eta + rlogis(n), c(-Inf, -1, 0, 1, Inf)), ordered = TRUE)
  d$ycat <- factor(sample(c("p", "q", "r"), n, TRUE))
  d$yz <- ifelse(runif(n) < 0.3, 0, d$yc)
  d$ybeta <- plogis(eta / 2 + rnorm(n, 0, 0.3))
  d$y2 <- rnorm(n, 0.4 * d$z)
  d$xm <- d$x + rnorm(n, 0, 0.2)
  d$xmi <- d$x
  d$xmi[c(5, 17, 40)] <- NA
  d
}
cases <- list(
  list("gauss", bf(y ~ x + z), gaussian()),
  list("inter", bf(y ~ x * f), gaussian()),
  list("smooth", bf(y ~ s(x) + z), gaussian()),
  list("gp", bf(y ~ gp(x, k = 8) + z), gaussian()),
  list("mo", bf(y ~ mo(o) + x), gaussian()),
  list("me", bf(y ~ me(xm, sdx) + z), gaussian()),
  list("mi", bf(y ~ mi(xmi) + z) + bf(xmi | mi() ~ z), gaussian()),
  list("nl", bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE),
       gaussian()),
  list("sigma", bf(y ~ x, sigma ~ z), gaussian()),
  list("ri", bf(y ~ x + (1 | g)), gaussian()),
  list("rs", bf(y ~ x + (1 + x | g)), gaussian()),
  list("crossed", bf(y ~ x + (1 | g) + (1 | h)), gaussian()),
  list("nested", bf(y ~ x + (1 | h / g)), gaussian()),
  list("grby", bf(y ~ x + (1 | gr(g, by = fg))), gaussian()),
  list("mm", bf(y ~ x + (1 | mm(g1, g2))), gaussian()),
  list("pois", bf(yc ~ x + (1 | g)), poisson()),
  list("binom", bf(yb | trials(tr) ~ x + (1 | g)), binomial()),
  list("cum", bf(yo ~ x + z), cumulative()),
  list("categ", bf(ycat ~ x), categorical()),
  list("zip", bf(yz ~ x, zi ~ z), zero_inflated_poisson()),
  list("hurdle", bf(yz ~ x), hurdle_poisson()),
  list("student", bf(y ~ x + z), student()),
  list("beta", bf(ybeta ~ x), Beta()),
  list("lognormal", bf(yp ~ x + (1 | g)), lognormal()),
  list("mv", bf(y ~ x) + bf(y2 ~ z), gaussian()),
  list("offset", bf(yc ~ x + offset(z / 10)), poisson()),
  list("fs", bf(y ~ x + s(z, g, bs = "fs", k = 4)), gaussian()),
  list("sigma_re", bf(y ~ x + (1 | g), sigma ~ (1 | g)), gaussian()),
  list("t2", bf(y ~ t2(x, z)), gaussian()),
  list("mmby", bf(y ~ x + (1 | mm(g1, g2, by = cbind(fg, fg)))), gaussian())
)
safe <- function(expr) {
  tryCatch(withCallingHandlers(expr, warning = function(w) {
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR:", conditionMessage(e)))
}
out <- list()
for (i in seq_along(cases)) {
  cs <- cases[[i]]
  nm <- cs[[1]]
  d <- mk(i)
  fit <- safe(frm(cs[[2]], family = cs[[3]], data = d))
  if (is.character(fit)) {
    out[[nm]] <- list(fit = fit)
    cat(nm, fit, "\n")
    next
  }
  resp <- names(fit$spec$responses)[1L]
  nd <- d[1:6, ]
  if ("g" %in% names(nd)) {
    nd$g <- factor(c(as.character(nd$g[1:4]), "new", "new"),
                   levels = c(levels(d$g), "new"))
  }
  r <- list()
  r$ce <- safe(conditional_effects(fit, resolution = 5))
  r$ce_null <- safe(conditional_effects(fit, resolution = 5,
                                        re_formula = NULL))
  r$ce_boot <- safe(conditional_effects(fit, resolution = 4, band = "boot",
                                        boot = 8, seed = 1))
  r$ce_boot_null <- safe(conditional_effects(fit, resolution = 4,
                                             band = "boot", boot = 8,
                                             seed = 1, re_formula = NULL))
  r$fitted <- safe(fitted(fit))
  r$fitted_nd <- safe(fitted(fit, newdata = nd, allow_new_levels = TRUE))
  set.seed(1)
  r$predict <- safe(predict(fit, ndraws = 50))
  set.seed(1)
  r$predict_nd <- safe(predict(fit, newdata = nd, allow_new_levels = TRUE,
                               ndraws = 50))
  set.seed(1)
  r$predict_ol <- safe(predict(fit, newdata = nd, allow_new_levels = TRUE,
                               sample_new_levels = "old_levels", ndraws = 50))
  out[[nm]] <- r
  cat(nm, "done\n")
}
saveRDS(out, sprintf(
  "C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-log/noreg-%s.rds",
  arm))
cat("saved\n")
