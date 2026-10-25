# The release review's findings 2, 3 and 6, probed on one build: what
# frmtmb_covstruct_info("sar") says with frmtmb.spatial loaded; the
# prior rows of a sar() model and of an ar() model; whether the
# internal row of lagsar takes a prior; and what class = "lagsar" on an
# error SAR says. 5 x 5 rook lattice, seed 2 (frmtmb.spatial's
# test-sar.R sar_data()).
#   Rscript dev/rel075-review-probe.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.spatial)})
cat("lib:", dirname(find.package("frmtmb")), format(packageVersion("frmtmb")),
    "\n")
say <- function(lab, expr) {
  r <- tryCatch({ v <- expr; if (is.data.frame(v)) "a table" else
    paste(format(v), collapse = " ") },
    error = function(e) paste0("ERROR [", paste(class(e)[1:2],
                                                collapse = ","), "] ",
                               conditionMessage(e)))
  cat(lab, ":", substr(r, 1, 400), "\n")
}
say("csinfo sar", frmtmb_covstruct_info("sar"))
say("csinfo fcor", frmtmb_covstruct_info("fcor"))
set.seed(2)
k <- 5
g <- expand.grid(r = 1:k, c = 1:k)
A <- 1 * (as.matrix(dist(g, method = "manhattan")) == 1)
W <- A / rowSums(A)
d <- data.frame(x = rnorm(k * k))
d$y <- as.vector(solve(diag(k * k) - 0.5 * W, 1 + 0.5 * d$x + rnorm(k * k)))
dp <- default_prior(y ~ x + sar(W), data = d, data2 = list(W = W))
print(as.data.frame(dp)[, c("prior", "class", "coef", "lb", "ub")])
d2 <- data.frame(t = 1:30, y = rnorm(30))
dpa <- default_prior(y ~ ar(time = t), data = d2)
print(as.data.frame(dpa)[, c("prior", "class", "coef")])
say("theta thetaac_1 on lag", validate_prior(
  set_prior("normal(0, 0.5)", class = "theta", coef = "thetaac_1"),
  y ~ x + sar(W), data = d, data2 = list(W = W)))
say("lagsar and theta both", validate_prior(
  c(set_prior("normal(0.3, 0.2)", class = "lagsar"),
    set_prior("normal(0, 0.5)", class = "theta", coef = "thetaac_1")),
  y ~ x + sar(W), data = d, data2 = list(W = W)))
say("lagsar on error", validate_prior(
  set_prior("normal(0, 0.5)", class = "lagsar"),
  y ~ x + sar(W, type = "error"), data = d, data2 = list(W = W)))
say("errorsar on error", validate_prior(
  set_prior("normal(0, 0.5)", class = "errorsar"),
  y ~ x + sar(W, type = "error"), data = d, data2 = list(W = W)))
