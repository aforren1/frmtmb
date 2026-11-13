# Probe the two predict() defects of lane ppdraw (backlog 3269, 3277).
# Usage: Rscript dev/ppdraw-probe1.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (length(args) && args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")

cat("\n== A. zero_inflated_binomial at the zi boundary (3269)\n")
n <- 200
set.seed(7)
dz <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE))
dz$y <- rbinom(n, dz$N, plogis(0.3 + 0.5 * dz$x))
fz <- withCallingHandlers(
  frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(), data = dz),
  message = function(m) cat("  M:", conditionMessage(m), "\n"),
  warning = function(w) cat("  W:", conditionMessage(w), "\n"))
print(fz$estimates)
print(summary(fz)$fixed)
V <- vcov(fz, full = TRUE)
print(sqrt(diag(V)))
sdr <- frmtmb:::sdr_of(fz)
cat("se_lost:\n"); print(sdr$se_lost)
set.seed(1)
p <- predict(fz, ndraws = 400, summary = FALSE)
cat("predict mean", mean(p), " fitted mean", mean(fitted(fz)[, "Estimate"]),
    "\n")
set.seed(1)
p0 <- predict(fz, ndraws = 400, summary = FALSE, propagate_error = FALSE)
cat("predict(propagate_error = FALSE) mean", mean(p0), "\n")

cat("\n== B. mixture(cumulative, cumulative) predict NA share (3277)\n")
set.seed(7)
do <- data.frame(x = rnorm(400))
ko <- rbinom(400, 1, 0.3)
do$yo <- cut(do$x + ifelse(ko == 1, 2, 0) + rlogis(400),
             c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
f4 <- withCallingHandlers(
  frm(bf(yo ~ x), family = mixture(cumulative(), cumulative()), data = do),
  message = function(m) cat("  M:", conditionMessage(m), "\n"),
  warning = function(w) cat("  W:", conditionMessage(w), "\n"))
print(f4$estimates)
print(sqrt(diag(vcov(f4, full = TRUE))))
set.seed(13)
nd <- data.frame(x = c(-1, 0, 1))
rp <- withCallingHandlers(
  predict(f4, newdata = nd, summary = FALSE, ndraws = 4000),
  warning = function(w) cat("  W:", conditionMessage(w), "\n"))
cat("NA share", mean(is.na(as.numeric(rp))), "\n")
print(table(as.numeric(rp), useNA = "always"))
set.seed(13)
rp0 <- predict(f4, newdata = nd, summary = FALSE, ndraws = 4000,
               propagate_error = FALSE)
cat("NA share, propagate_error = FALSE", mean(is.na(as.numeric(rp0))), "\n")
set.seed(13)
rpi <- predict(f4, summary = FALSE, ndraws = 400)
cat("NA share in-sample", mean(is.na(as.numeric(rpi))), "\n")
print(fitted(f4, newdata = nd))
