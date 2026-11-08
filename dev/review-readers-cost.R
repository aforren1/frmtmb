# Review of lane readers: the cost of the integral at a new level, and what
# marginaleffects pays for it. REV_ARM = before | after | mfxse (the mfxse
# lane's library ahead of rellib-r18, read-only). Each timing is the minimum
# of 5 rounds. Data seed 4101 as dev/readers-item1.R.
arm <- Sys.getenv("REV_ARM", "after")
if (arm == "mfxse") {
  .libPaths(c("C:/Users/adf44/source/r/wt-mfxse-lib",
              "C:/Users/adf44/source/r/rellib-r18",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressPackageStartupMessages(library(frmtmb))
  cat("arm: mfxse lib:", dirname(find.package("frmtmb")), "\n")
  ns <- asNamespace("frmtmb")
} else {
  source("dev/review-readers-setup.R")
}
suppressPackageStartupMessages(library(marginaleffects))
cat("marginaleffects", format(packageVersion("marginaleffects")), "\n")
tmin <- function(f, k = 5) {
  min(replicate(k, system.time(f())[["elapsed"]]))
}
set.seed(4101)
J <- 40; n <- 400
g <- factor(sample(seq_len(J), n, TRUE))
u <- rnorm(J, 0, 1)
x <- runif(n, -1, 1)
dd <- data.frame(x = x, g = g)
dd$y <- rpois(n, exp(0.5 + 0.3 * x + u[g]))
fit <- frm(bf(y ~ x + (1 | g)), family = poisson(), data = dd)
nd3 <- data.frame(x = c(-0.5, 0, 0.8), g = "new")
nd100 <- data.frame(x = seq(-1, 1, length.out = 100), g = "new")
cat("fitted, 3 rows new level:  ",
    tmin(function() fitted(fit, newdata = nd3, allow_new_levels = TRUE)), "s\n")
cat("fitted, 100 rows new level:",
    tmin(function() fitted(fit, newdata = nd100, allow_new_levels = TRUE)), "s\n")
cat("frm_lp_basis, 100 rows (a closed form's one call):",
    tmin(function() frm_lp_basis(fit, newdata = nd100, allow_new_levels = TRUE)), "s\n")
# count get_predict() calls of one avg_predictions(): a counting wrapper
# registered over frmtmb's method
cnt_env <- new.env()
cnt_env$n <- 0L
orig <- get("get_predict.frmtmb_fit", ns)
registerS3method("get_predict", "frmtmb_fit", function(model, newdata, ...) {
  cnt_env$n <- cnt_env$n + 1L
  orig(model, newdata, ...)
}, envir = asNamespace("marginaleffects"))
ap_run <- function(label, nd) {
  cnt_env$n <- 0L
  t <- system.time(ap <- suppressWarnings(
    avg_predictions(fit, newdata = nd, allow_new_levels = TRUE)))[["elapsed"]]
  a <- as.data.frame(ap)
  cat(sprintf("%-40s %7.2f s %4d get_predict() calls; estimate %.5f se %.5f\n",
              label, t, cnt_env$n, a$estimate, a$std.error))
}
ap_run("avg_predictions, 100 new-level rows", nd100)
ap_run("avg_predictions, 400 seen rows, flag on", dd)
nd_s <- dd[1:100, ]
nd_s$g <- as.character(nd_s$g)
nd_s$g[1] <- "new"
ap_run("avg_predictions, 99 seen + 1 new row", nd_s)
