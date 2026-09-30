# Reviewer: frm_sample() draws of models without an equation must be
# identical() between arms; with an equation, sigma1 must equal sigma2
# on every draw. Usage: Rscript dev/formula2-rev-sample.R before|after
arm <- commandArgs(TRUE)[1]
base <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "after") c("C:/Users/adf44/source/r/wt-formula2-lib",
                                base) else base)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb.sample from", find.package("frmtmb.sample"), "\n")
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                f2 = factor(sample(c("p", "q"), n, TRUE)))
k <- rbinom(n, 1, 0.4)
d$y <- ifelse(k == 1, rnorm(n, 3 + 0.5 * d$x, 1), rnorm(n, -1, 1))
fam <- mixture(gaussian(), gaussian())
q <- function(expr) suppressWarnings(suppressMessages(expr))
out <- list()
fits <- list(
  gs = q(frm(bf(y ~ x, sigma ~ f2), data = d)),
  mix = q(frm(bf(y ~ x), family = fam, data = d))
)
for (nm in names(fits)) {
  s <- q(frm_sample(fits[[nm]], chains = 1, iter = 300, warmup = 150,
                    refresh = 0, seed = 3))
  dr <- as.matrix(posterior::as_draws_matrix(s))
  out[[nm]] <- list(draws = dr, prior = q(prior_summary(s)))
  cat(nm, "draws", dim(dr), "\n")
}
if (arm == "after") {
  fe <- q(frm(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d))
  s <- q(frm_sample(fe, chains = 2, iter = 300, warmup = 150, refresh = 0,
                    seed = 3))
  dr <- as.matrix(posterior::as_draws_matrix(s))
  cat("equated columns:", colnames(dr), "\n")
  cat("sigma1 identical to sigma2 on every draw:",
      identical(dr[, "sigma1"], dr[, "sigma2"]), " n draws", nrow(dr), "\n")
  cat("variables():", variables(s), "\n")
  print(q(prior_summary(s)))
  print(q(summary(s)))
  cat("posterior_summary rows:", rownames(q(posterior_summary(s))), "\n")
  fx <- q(fixef(s)); print(fx)
  ft <- q(fitted(s, dpar = "sigma1", newdata = d[1:3, ]))
  ft2 <- q(fitted(s, dpar = "sigma2", newdata = d[1:3, ]))
  cat("fitted sigma1 == sigma2 draws summary:", identical(ft, ft2), "\n")
  pp <- q(posterior_predict(s, ndraws = 20))
  cat("posterior_predict dim", dim(pp), "\n")
}
saveRDS(out, sprintf(
  "C:/Users/adf44/source/r/frmtmb-wt-formula2/dev/formula2-rev-sample-%s.rds",
  arm))
cat("DONE\n")
