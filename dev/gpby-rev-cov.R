# Reviewer: coverage of frmtmb.spline's simultaneous band over an exact
# gp(), the lane's design (dev/gpby-cov-band.R: 60 points on [0, 6], a
# 25-point grid on [4, 9], f one draw of the model's own GP with sd 1 and
# length scale 1, y = 0.5 + f + N(0, 0.2^2), nsim = 4000 with the seed
# as the simulation seed), in two modes:
#   plugin: the lane's fit, y ~ gp(x) by maximum likelihood (reproduces
#           the lane's per-seed rows);
#   oracle: the same fit with the GP sd, its length scale and sigma held
#           at the truth by priors of sd 1e-4 (internal theta) and 1e-5
#           (sigma), so the band is judged on its own calibration and not
#           on the plug-in error of the hyperparameters.
# Usage: Rscript dev/gpby-rev-cov.R <arm> <mode> <from> <to>
args <- commandArgs(TRUE)
arm <- args[1]; mode <- args[2]
s_from <- as.integer(args[3]); s_to <- as.integer(args[4])
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
dir.create(file.path(wt, "dev/gpby-rev-cov"), showWarnings = FALSE)
out <- file.path(wt, sprintf("dev/gpby-rev-cov/%s-%s-%d-%d.tsv", arm, mode,
                             s_from, s_to))
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spline"), "\n",
    file = paste0(out, ".lib"))
grid <- seq(4, 9, length.out = 25)
pr <- set_prior("normal(0, 1e-4)", class = "theta", coef = "theta_1") +
  set_prior("normal(0, 1e-4)", class = "theta", coef = "theta_2") +
  set_prior("normal(0.2, 1e-5)", class = "sigma")
cat("seed\tok\tsim_cover\tpt_cover\tcrit\tse_max\tcf_rel\n", file = out)
for (s in s_from:s_to) {
  set.seed(s)
  x <- sort(stats::runif(60, 0, 6))
  allx <- c(x, grid)
  K <- exp(-outer(allx, allx, "-")^2 / 2)
  f <- drop(t(chol(K + diag(1e-8, length(allx)))) %*%
              stats::rnorm(length(allx)))
  d <- data.frame(x = x, y = 0.5 + f[seq_along(x)] +
                    stats::rnorm(60, 0, 0.2))
  truth <- 0.5 + f[length(x) + seq_along(grid)]
  r <- tryCatch({
    fit <- suppressWarnings(if (mode == "oracle") {
      frm(bf(y ~ gp(x)), data = d, prior = pr)
    } else frm(bf(y ~ gp(x)), data = d))
    cv <- suppressWarnings(frm_curve(fit, newdata = data.frame(x = grid),
                                     nsim = 4000, seed = s))
    # the closed form at the oracle's hyperparameters: the posterior
    # covariance of 0.5 + f(grid) given y with a flat intercept, with
    # frmtmb's 1e-6 nugget on the kernel
    cf <- NA_real_
    if (mode == "oracle") {
      ux <- sort(unique(x))
      kf <- function(a, b) exp(-outer(a, b, "-")^2 / 2) +
        1e-6 * (outer(a, b, "-") == 0)
      Z <- outer(x, ux, "==") * 1
      C <- Z %*% kf(ux, ux) %*% t(Z) + diag(0.04, 60)
      Ci <- solve(C)
      Kgx <- kf(grid, ux) %*% t(Z)
      R <- 1 - Kgx %*% Ci %*% rep(1, 60)
      S <- kf(grid, grid) - Kgx %*% Ci %*% t(Kgx) +
        R %*% t(R) / sum(Ci)
      Sg <- attr(cv, "Sigma")
      cf <- max(abs(Sg - S)) / max(abs(S))
    }
    c(ok = 1,
      sim = as.numeric(all(truth >= cv$.lower_sim &
                             truth <= cv$.upper_sim)),
      pt = mean(truth >= cv$.lower_ci & truth <= cv$.upper_ci),
      crit = cv$.crit_sim[1], se = max(cv$.se), cf = cf)
  }, error = function(e) c(ok = 0, sim = NA, pt = NA, crit = NA, se = NA,
                           cf = NA))
  cat(sprintf("%d\t%d\t%s\t%s\t%s\t%s\t%s\n", s, r[["ok"]], r[["sim"]],
              r[["pt"]], r[["crit"]], r[["se"]], r[["cf"]]),
      file = out, append = TRUE)
}
cat("DONE\n", file = paste0(out, ".lib"), append = TRUE)
