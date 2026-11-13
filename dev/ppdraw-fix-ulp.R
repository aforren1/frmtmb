# Fix round: do the new assertions of test-predict-draw-edge.R hold on
# every rounding? The covariate (or response, for the gaussian fit) is
# perturbed by x * (1 + 2e-15 * k), k = 0..20, and each fit is asked what
# the test asks. Usage: Rscript dev/ppdraw-fix-ulp.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
said <- function(expr) {
  m <- character(0)
  v <- withCallingHandlers(expr, frmtmb_draw_held = function(c) {
    m <<- c(m, conditionMessage(c))
    invokeRestart("muffleMessage")
  }, warning = function(w) invokeRestart("muffleWarning"))
  list(v = v, m = paste(m, collapse = " "))
}
zf <- function(fit) {
  set.seed(1)
  r <- said(predict(fit, ndraws = 400, summary = FALSE))
  m <- rowMeans(r$v)
  list(z = (mean(m) - mean(fitted(fit)[, "Estimate"])) /
         (stats::sd(m) / sqrt(400)), zero = mean(r$v == 0), m = r$m)
}
for (k in 0:20) {
  e <- 1 + 2e-15 * k
  set.seed(7)
  dz <- data.frame(x = stats::rnorm(200), N = sample(5:15, 200, TRUE))
  dz$y <- stats::rbinom(200, dz$N, stats::plogis(0.3 + 0.5 * dz$x))
  dz$x <- dz$x * e
  f1 <- suppressWarnings(frm(bf(y | trials(N) ~ x, zi ~ x),
                             family = zero_inflated_binomial(), data = dz))
  a1 <- zf(f1)
  set.seed(12)
  dh <- data.frame(x = stats::rnorm(300))
  dh$y <- stats::rpois(300, exp(2.5 + 0.2 * dh$x))
  dh$x <- dh$x * e
  f4 <- suppressWarnings(frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(),
                             data = dh))
  a4 <- zf(f4)
  set.seed(3)
  dp <- data.frame(x = stats::rnorm(200))
  dp$y <- stats::rpois(200, exp(1 + 0.3 * dp$x))
  dp$x <- dp$x * e
  f3 <- suppressWarnings(frm(bf(y ~ x), family = zero_inflated_poisson(),
                             data = dp))
  a3 <- zf(f3)
  set.seed(5)
  dg <- data.frame(g = factor(rep(1:6, each = 20)), x = stats::rnorm(120))
  dg$y <- (1 + 0.5 * dg$x + stats::rnorm(6, 0, 0.1)[dg$g] +
             stats::rnorm(120)) * e
  fg <- suppressWarnings(frm(bf(y ~ x + (1 | g)), data = dg))
  set.seed(1)
  pg <- said(predict(fg, newdata = data.frame(x = 0, g = factor("new")),
                     ndraws = 4000, summary = FALSE,
                     allow_new_levels = TRUE))
  set.seed(1)
  p0 <- predict(fg, newdata = data.frame(x = 0, g = factor("new")),
                ndraws = 4000, summary = FALSE, allow_new_levels = TRUE,
                propagate_error = FALSE)
  cat(sprintf(paste("k=%2d A1 z %6.2f narrowed %s | A4 z %6.2f narrowed %s |",
                    "zip3 zi %.3f z %6.2f zero %.3f narrowed %s flat %s |",
                    "g6 sd ratio %.3f flat %s\n"), k,
              a1$z, grepl("_x` are narrowed", a1$m), a4$z,
              grepl("_x` are narrowed", a4$m), f3$estimates$betad[1], a3$z,
              a3$zero, grepl("`zi_(Intercept)` are narrowed", a3$m,
                             fixed = TRUE),
              grepl("flat toward", a3$m), stats::sd(pg$v) / stats::sd(p0),
              grepl("flat toward", pg$m)))
}
