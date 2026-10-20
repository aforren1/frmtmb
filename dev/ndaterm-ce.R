# Lane ndaterm: conditional_effects(method = "predict") on a wiener
# dec() model. What value does the grid hold the boundary at, and does a
# factor dec() reach the simulator coerced as the fit coerced it?
# Seed 7. Usage: Rscript dev/ndaterm-ce.R <arm>; log
# dev/ndaterm-out/ce-<arm>.txt
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "r9"
base <- "C:/Users/adf44/source/r"
libs <- switch(arm,
  mine = c(file.path(base, "wt-ndaterm-lib"), file.path(base, "rellib-r9")),
  file.path(base, paste0("rellib-", arm)))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
set.seed(7)
dat <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.3)
dat$x <- rnorm(300)
dat$side <- factor(ifelse(dat$upper == 1, "up", "down"),
                   levels = c("up", "down"))
cat("mean(upper)", mean(dat$upper), "\n")
fw <- frm(bf(rt | dec(upper) ~ x, bias = 0.3), family = wiener(),
          data = dat)
ff <- frm(bf(rt | dec(side) ~ x, bias = 0.3), family = wiener(),
          data = dat)
cat("loglik numeric dec", logLik(fw), " factor dec", logLik(ff),
    " (levels up, down: the coercion codes 'up' as 0, the lower ",
    "boundary)\n")
med <- function(fit, nd) {
  s <- simulate(fit, nsim = 2000, newdata = nd, seed = 7)
  apply(as.matrix(s), 1, median)
}
nd0 <- data.frame(x = 0, upper = c(1, 0))
cat("simulate median at x = 0, upper = 1, 0:",
    format(med(fw, nd0), digits = 4), "\n")
one <- function(fit, cond = NULL) {
  set.seed(7)
  a <- list(fit, "x", method = "predict", ndraws = 2001)
  if (!is.null(cond)) a$conditions <- cond
  ce <- do.call(conditional_effects, a)
  d <- ce[[1]]
  i <- which.min(abs(d$x))
  held <- intersect(c("upper", "side"), names(d))
  sprintf("x = %.3f, held %s = %s, estimate %.4f", d$x[i], held,
          format(d[[held]][i]), d$estimate__[i])
}
cat("numeric dec, default grid:",
    tryCatch(one(fw), error = conditionMessage), "\n")
cat("numeric dec, upper = 1:",
    tryCatch(one(fw, data.frame(upper = 1)), error = conditionMessage), "\n")
cat("numeric dec, upper = 0:",
    tryCatch(one(fw, data.frame(upper = 0)), error = conditionMessage), "\n")
cat("factor dec, default grid:",
    tryCatch(one(ff), error = conditionMessage), "\n")
cat("factor dec, side = up:",
    tryCatch(one(ff, data.frame(side = factor("up", c("up", "down")))),
             error = conditionMessage), "\n")
cat("factor dec, side = down:",
    tryCatch(one(ff, data.frame(side = factor("down", c("up", "down")))),
             error = conditionMessage), "\n")
ndf <- data.frame(x = 0, side = factor(c("up", "down"), c("up", "down")))
cat("factor dec simulate median at side = up, down:",
    format(med(ff, ndf), digits = 4), "\n")
cat("factor dec fitted at side = up, down:",
    format(fitted(ff, newdata = ndf)[, "Estimate"], digits = 4), "\n")
cat("numeric dec fitted at upper = 0, 1:",
    format(fitted(fw, newdata = data.frame(x = 0, upper = c(0, 1)))[,
      "Estimate"], digits = 4), "\n")

cat("\n== custom family with vint() trials, conditional_effects(predict)\n")
cbin <- custom_family(
  "cbin", dpars = "mu", links = list(mu = "logit"),
  lpdf = function(y, dpars, aterms) {
    RTMB::dbinom(y, aterms[["vint1"]], dpars[["mu"]], log = TRUE)
  },
  required_aterms = "vint1", type = "discrete",
  post = list(mean_fn = function(dpars, aterms) {
    dpars[["mu"]] * aterms[["vint1"]]
  }),
  sim = function(dpars, aterms, n) {
    stats::rbinom(n, aterms[["vint1"]], dpars[["mu"]])
  })
set.seed(7)
db <- data.frame(x = rnorm(200), N = sample(5:15, 200, TRUE))
db$y <- rbinom(200, db$N, plogis(0.3 + 0.5 * db$x))
fc <- frm(bf(y | vint(N) ~ x), family = cbin, data = db)
cat("mean(N)", mean(db$N), " round(median(N))", round(median(db$N)), "\n")
r <- tryCatch(withCallingHandlers({
  set.seed(7)
  d <- conditional_effects(fc, "x", method = "predict")[[1]]
  sprintf("held N = %s, estimate range %s", format(d$N[1]),
          paste(format(range(d$estimate__)), collapse = " to "))
}, warning = function(w) {
  cat("  warning:", conditionMessage(w), "\n")
  invokeRestart("muffleWarning")
}), error = function(e) paste("ERROR:", conditionMessage(e)))
cat("vint ce predict, default grid:", r, "\n")
