# Coverage of frmtmb.spline's simultaneous band over an exact gp(), on a
# grid that runs from inside the data to past its edge.
#
# Per seed: the true latent field f is ONE draw of the model's own GP
# (sd 1, length scale 1) at the 60 observed positions and the 25 grid
# positions jointly; y = 0.5 + f(x) + N(0, 0.2^2). The band covers when
# 0.5 + f(grid) is inside it at every grid point. Pointwise coverage is
# recorded beside it.
#
# Usage: Rscript dev/gpby-cov-band.R <arm: base|lane> <seed_from> <seed_to>
# Writes one line per seed to dev/gpby-cov/<arm>-<from>-<to>.tsv.
args <- commandArgs(TRUE)
arm <- args[1]
s_from <- as.integer(args[2])
s_to <- as.integer(args[3])
LIB <- "C:/Users/adf44/source/r/wt-gpby-lib"
.libPaths(c(if (arm == "lane") LIB, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
dir.create("dev/gpby-cov", showWarnings = FALSE)
out <- sprintf("dev/gpby-cov/%s-%d-%d.tsv", arm, s_from, s_to)
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spline"), "\n",
    file = paste0(out, ".lib"))
grid <- seq(4, 9, length.out = 25)
cat("seed\tok\tsim_cover\tpt_cover\tcrit\tse_max\n", file = out)
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
    fit <- suppressWarnings(frm(bf(y ~ gp(x)), data = d))
    cv <- suppressWarnings(frm_curve(fit, newdata = data.frame(x = grid),
                                     nsim = 4000, seed = s))
    c(ok = 1,
      sim = as.numeric(all(truth >= cv$.lower_sim &
                             truth <= cv$.upper_sim)),
      pt = mean(truth >= cv$.lower_ci & truth <= cv$.upper_ci),
      crit = cv$.crit_sim[1], se = max(cv$.se))
  }, error = function(e) c(ok = 0, sim = NA, pt = NA, crit = NA, se = NA))
  cat(sprintf("%d\t%d\t%s\t%s\t%s\t%s\n", s, r[["ok"]], r[["sim"]],
              r[["pt"]], r[["crit"]], r[["se"]]),
      file = out, append = TRUE)
}
cat("DONE\n", file = paste0(out, ".lib"), append = TRUE)
