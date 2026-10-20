# Reviewer of lane ndaterm, fix round 1: frmtmb.sample's
# posterior_predict() of a gaussian mixture on newdata and under
# re_formula = NA, against posterior_epred(). Usage: <arm>. Seed 7.
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "mine"
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ndaterm/dev/stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, find.package("frmtmb.sample"), "\n")
set.seed(7)
n <- 300
d <- data.frame(x = rnorm(n))
k <- rbinom(n, 1, plogis(-1.4 + 0.8 * d$x))
d$y <- ifelse(k == 1, 6, 0) + d$x + rnorm(n, 0, 0.5)
f1 <- frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = d)
dr <- suppressWarnings(suppressMessages(frm_sample(
  f1, chains = 1, iter = 600, warmup = 300, refresh = 0, seed = 7)))
nd <- data.frame(x = c(-1, 0, 1))
set.seed(1)
p <- posterior_predict(dr, newdata = nd)
e <- posterior_epred(dr, newdata = nd)
se <- apply(p - e, 2, sd) / sqrt(nrow(p))
cat("newdata: predict means", format(colMeans(p), digits = 4),
    "\n         epred means  ", format(colMeans(e), digits = 4),
    "\n         6 se         ", format(6 * se, digits = 4), "\n")
set.seed(2)
p0 <- posterior_predict(dr, re_formula = NA)
e0 <- posterior_epred(dr)
dd <- colMeans(p0) - colMeans(e0)
cat("in sample re_formula = NA: mean |predict - epred| / (sd/sqrt(n))",
    format(mean(abs(dd) / (apply(p0 - e0, 2, sd) / sqrt(nrow(p0)))),
           digits = 4), "(about 0.8 when unbiased)\n")
cat("done\n")
