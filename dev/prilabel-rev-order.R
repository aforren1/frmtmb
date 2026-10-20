# Review of lane prilabel: variables() order beside smooths and gp(),
# and ranef() names with two blocks on one factor, on both arms.
#   arg 1: "lane" or "base". Data seed 61.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
suppressPackageStartupMessages(library(frmtmb))
sink(file.path(wt, paste0("dev/prilabel-rev-log/order-", arm, ".txt")),
     split = TRUE)
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(61)
n <- 240
d <- data.frame(x = runif(n), g = factor(sample(1:10, n, TRUE)),
                h = factor(sample(1:6, n, TRUE)))
d$y <- sin(4 * d$x) + rnorm(10, 0, 0.6)[d$g] + rnorm(10, 0, 0.4)[d$g] * d$x +
  rnorm(6, 0, 0.4)[d$h] + rnorm(n, 0, 0.3)
for (f in list(y ~ s(x) + (1 + x | g),
               y ~ (1 + x | g) + s(x),
               y ~ (1 + x | h) + (1 + x | g),
               y ~ (1 | g) + (0 + x | g))) {
  ft <- frm(bf(f), data = d)
  cat("\n--", deparse1(f), "\n")
  v <- variables(ft)
  print(v[!grepl("^b_|^r_|^s_|^zs_|^lprior|^lp__", v)])
  cat("names(ranef):", names(ranef(ft)), "\n")
}
sink()
if (arm == "lane") {
  .libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
  suppressPackageStartupMessages(library(brms))
  sink(file.path(wt, "dev/prilabel-rev-log/order-brms.txt"), split = TRUE)
  b <- brm(y ~ (1 | h) + s(x) + (1 + x | g) + (0 + x | h), data = d,
           chains = 1, iter = 100, warmup = 50, seed = 1, refresh = 0)
  v <- variables(b)
  print(v[!grepl("^b_|^r_|^s_|^zs_|^lprior|^lp__|^bs_", v)])
  cat("names(ranef):", names(ranef(b)), "\n")
  print(lapply(ranef(b), function(a) dimnames(a)[[3]]))
  sink()
}
