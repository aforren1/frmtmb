# Reviewer of lane priortab: a default's bound kept for resp = "Ya" vs "Y_a"
# (student nu above 1 on the sampling route). Lane library only.
.libPaths(c("C:/Users/adf44/source/r/wt-priortab-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1)
n <- 120
d <- data.frame(z = rnorm(n))
d$Y_a <- 1 + d$z * 0.3 + rt(n, 5) * 0.3
d$Y_b <- d$Y_a + rnorm(n)
f <- bf(mvbind(Y_a, Y_b) ~ z) + set_rescor(FALSE)
fit <- frm(f, data = d, family = student(), dry_run = "objective")
for (r in c("Y_a", "Ya")) {
  pl <- set_prior("gamma(2, 0.1)", class = "nu", resp = r)
  rp <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, pl))
  cat("\nresp =", r, "\n")
  cat("ri$lower:\n"); print(rp$ri$lower)
  x <- as.data.frame(rp$effective)
  print(x[x$class == "nu", intersect(c("prior", "class", "resp", "lb"),
                                     names(x))])
}
