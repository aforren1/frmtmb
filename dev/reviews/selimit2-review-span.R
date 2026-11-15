# Review of lane selimit2, item 5433: is "inside the fitted span" a proof
# that a row's mean stays put over a limit's level range? zi ~ x + z: rows
# A (x = 0, z in [0, 1]) fix zi's intercept and z slope; rows B (x in
# [1, 2], z = 0) have no zeros, so zi_x runs to -Inf. A new row at x = 1,
# z = 1 lies inside each column's observed range and inside the span of
# zi_x's move over the fitted rows ([0, 2]); the likelihood is level for
# zi_x down to about -9, where that row's zi is 0.4.
#   Rscript selimit2-review-span.R <base|lane>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "lib", dirname(find.package("frmtmb")), "\n")
set.seed(31)
nA <- 200
nB <- 100
zA <- stats::runif(nA)
ziA <- stats::plogis(-9 + 18 * zA)
yA <- ifelse(stats::runif(nA) < ziA, 0L, stats::rpois(nA, 4))
xB <- stats::runif(nB, 1, 2)
yB <- 1L + stats::rpois(nB, 3)
d <- data.frame(x = c(rep(0, nA), xB), z = c(zA, rep(0, nB)), y = c(yA, yB))
ws <- character(0)
fit <- withCallingHandlers(
  frm(bf(y ~ 1, zi ~ x + z), data = d, family = zero_inflated_poisson()),
  warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    cat("M:", substr(conditionMessage(m), 1, 250), "\n")
    invokeRestart("muffleMessage")
  })
for (w in ws) cat("W:", substr(w, 1, 250), "\n")
p <- fit$opt$par
cat("code", fit$opt$convergence, "par", paste(names(p), format(p, digits = 6)),
    "\n")
cat("se_lost:", paste(names(sdr_of <- frmtmb:::sdr_of(fit)$se_lost),
                      sdr_of), "\n")
nd <- data.frame(x = c(1, 1, 1.5, 0, 3), z = c(1, 0.8, 0.9, 0.5, 0))
fw <- character(0)
fp <- withCallingHandlers(stats::fitted(fit, newdata = nd),
                          warning = function(w) {
                            fw <<- c(fw, conditionMessage(w))
                            invokeRestart("muffleWarning")
                          })
for (w in fw) cat("fitted W:", substr(w, 1, 250), "\n")
print(cbind(nd, fp[, c("Estimate", "Est.Error")]))
# the level range of zi_x, the rest held, and each new row's mean there
obj <- fit$obj
f0 <- obj$fn(p)
jx <- which(names(p) == "betad")[2]
jz <- which(names(p) == "betad")[3]
ji <- which(names(p) == "betad")[1]
mu <- exp(p[which(names(p) == "beta")[1]])
for (bx in c(p[jx], -30, -20, -15, -12, -10, -9.5, -9, -8, -6)) {
  q <- p
  q[jx] <- bx
  dl <- obj$fn(q) - f0
  zi <- stats::plogis(q[ji] + bx * nd$x + q[jz] * nd$z)
  cat(sprintf("zi_x %8.3f | objective change %.3g | means %s\n", bx, dl,
              paste(format((1 - zi) * mu, digits = 5), collapse = " ")))
}
