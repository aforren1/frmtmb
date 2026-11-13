# Fix round: the acceptance cases of the coordinator's decision 1d, on one
# build. predict() against fitted() (z in Monte Carlo standard errors of
# the draws' mean) and P(y = 0), quantile widths where brms is the
# reference. Usage: Rscript dev/ppdraw-fix-cases.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
q <- function(expr) suppressMessages(suppressWarnings(expr))
said <- function(expr) {
  m <- character(0)
  v <- withCallingHandlers(expr, frmtmb_draw_held = function(c) {
    m <<- c(m, conditionMessage(c))
    invokeRestart("muffleMessage")
  }, warning = function(w) invokeRestart("muffleWarning"))
  list(v = v, m = m)
}
pz <- function(fit, ndraws = 400L, seed = 1, newdata = NULL) {
  set.seed(seed)
  r <- said(predict(fit, newdata = newdata, ndraws = ndraws,
                    summary = FALSE))
  p <- r$v
  m <- rowMeans(p)
  fm <- mean(q(fitted(fit, newdata = newdata))[, "Estimate"])
  list(z = (mean(m) - fm) / (stats::sd(m) / sqrt(ndraws)), pm = mean(m),
       fm = fm, p0 = mean(p == 0), m = r$m, p = p)
}
widths <- function(p) {
  c(w90 = mean(apply(p, 2, stats::quantile, 0.95) -
                  apply(p, 2, stats::quantile, 0.05)),
    w98 = mean(apply(p, 2, stats::quantile, 0.99) -
                  apply(p, 2, stats::quantile, 0.01)))
}
cat("\n== 3269 zib, 4000 draws (brms: mean 5.7022, P(y=0) 0.0097,",
    "w90 4.9752, w98 7.0005)\n")
set.seed(7)
dz <- data.frame(x = stats::rnorm(200), N = sample(5:15, 200, TRUE))
dz$y <- stats::rbinom(200, dz$N, stats::plogis(0.3 + 0.5 * dz$x))
fz <- q(frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
            data = dz))
r <- pz(fz, 4000L)
cat(sprintf("mean %.4f fitted %.4f z %.2f P(y=0) %.4f w90 %.4f w98 %.4f\n",
            r$pm, r$fm, r$z, r$p0, widths(r$p)[1], widths(r$p)[2]))
cat("message:", r$m, "\n")

cat("\n== A1 zib zi ~ x; A4 hurdle_poisson hu ~ x (400 draws)\n")
f1 <- q(frm(bf(y | trials(N) ~ x, zi ~ x), family = zero_inflated_binomial(),
            data = dz))
r <- pz(f1)
cat(sprintf("A1 mean %.4f fitted %.4f z %.2f\n", r$pm, r$fm, r$z))
cat("message:", r$m, "\n")
set.seed(12)
dh <- data.frame(x = stats::rnorm(300))
dh$y <- stats::rpois(300, exp(2.5 + 0.2 * dh$x))
f4 <- q(frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(), data = dh))
r <- pz(f4)
cat(sprintf("A4 mean %.4f fitted %.4f z %.2f\n", r$pm, r$fm, r$z))
cat("message:", r$m, "\n")

cat("\n== ZIP null data, seeds 1..30, 200 rows (brms seed 3: mean 3.0632,",
    "P(y=0) 0.0802)\n")
for (seed in 1:30) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.3 * d$x))
  fit <- q(frm(bf(y ~ x), family = zero_inflated_poisson(), data = d))
  r <- pz(fit, if (seed == 3) 4000L else 400L)
  cat(sprintf("seed %2d zi %8.3f mean %.4f fitted %.4f z %7.2f P(y=0) %.4f %s\n",
              seed, fit$estimates$betad[1], r$pm, r$fm, r$z, r$p0,
              if (length(r$m)) substr(r$m, 13, 70) else "-"))
}

cat("\n== gaussian (1 | g), 6 groups x 20, seed 5: new level, 4000 draws\n")
set.seed(5)
d <- data.frame(g = factor(rep(1:6, each = 20)), x = stats::rnorm(120))
d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, 0.1)[d$g] + stats::rnorm(120)
fit <- q(frm(bf(y ~ x + (1 | g)), data = d))
nd <- data.frame(x = 0, g = factor("new"))
for (pe in c(TRUE, FALSE)) {
  set.seed(1)
  r <- said(predict(fit, newdata = nd, ndraws = 4000, summary = FALSE,
                    allow_new_levels = TRUE, propagate_error = pe))
  cat(sprintf("propagate_error %-5s sd %.4g q01 %.4g q99 %.4g %s\n", pe,
              stats::sd(r$v), stats::quantile(r$v, 0.01),
              stats::quantile(r$v, 0.99),
              if (length(r$m)) substr(r$m, 1, 120) else ""))
}

cat("\n== se(s, sigma = TRUE), rows s = 0.05 and 5, 4000 draws",
    "(brms column sd 0.1077 5.0106)\n")
set.seed(7)
dg <- data.frame(x = stats::rnorm(200), s = stats::runif(200, 0.2, 1))
dg$y <- 1 + dg$x + stats::rnorm(200, 0, dg$s)
fg <- q(frm(bf(y | se(s, sigma = TRUE) ~ x), family = gaussian(), data = dg))
set.seed(1)
r <- said(predict(fg, newdata = data.frame(x = c(0, 0), s = c(0.05, 5)),
                  ndraws = 4000, summary = FALSE))
cat("column sd:", format(apply(r$v, 2, stats::sd), digits = 4), "\n")
