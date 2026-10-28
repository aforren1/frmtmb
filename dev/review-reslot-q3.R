# Reviewer short check: what frm_sample() says when handed frm()'s route
# table (default_prior()'s default), edited in one row as a brms user
# would. Usage: Rscript dev/review-reslot-q3.R before|after. Seed 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
set.seed(1)
G <- 10L; n <- 8L
d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n))
d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 0.7)[d$g] + rnorm(G * n)
f1 <- y ~ x + (1 + x | g)
p <- default_prior(f1, data = d)
cat("route attribute:", attr(p, "route"), "\n")
p$prior[p$class == "b" & p$coef == "x"] <- "normal(0, 1)"
cat("route attribute after the edit:", attr(p, "route"), "\n")
fitd <- suppressMessages(frm(f1, data = d, dry_run = "objective"))
withCallingHandlers(
  r <- frmtmb.sample:::sample_resolve_priors(fitd, p),
  message = function(m) { cat("MESSAGE:", conditionMessage(m)); invokeRestart("muffleMessage") },
  warning = function(w) { cat("WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
print(as.data.frame(r$effective)[, c("prior", "class", "coef", "group")],
      row.names = FALSE)
