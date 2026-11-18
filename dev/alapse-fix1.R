# Lane alapse, fix round 1: the review's cases on the fixed build.
#   Rscript dev/alapse-fix1.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/alapse-common.R")
show <- function(id, r) {
  f <- r$v
  if (inherits(f, "error")) { cat(id, "ERROR", conditionMessage(f), "\n"); return() }
  ce <- f$opt$ald_cert
  cat(sprintf("%-22s code %d ok %-5s global %-5s ll %.6f | %s | q %s | %s\n", id,
              f$opt$convergence, isTRUE(ce$ok), isTRUE(ce$global),
              as.numeric(logLik(f)),
              paste(unique(sub("^frmtmb_", "", r$wc)), collapse = ","),
              format(signif(ns$eval_dpars(f)[[1]]$quantile[1], 4)),
              ce$why %||% ce$local %||% ""))
  if (!is.null(ce$search)) cat("   search:", paste(names(ce$search), signif(unlist(ce$search), 5), collapse = " "), "\n")
}
set.seed(2000 + 109); d <- data.frame(x = rnorm(100)); d$y <- 1 + 0.5 * d$x + rnorm(100)
show("B109", cap(frm(y ~ x, family = asym_laplace(), data = d)))
for (s in c(22, 140)) {
  set.seed(2000 + s); d <- data.frame(x = rnorm(100)); d$y <- 1 + 0.5 * d$x + rnorm(100)
  show(paste0("B", s), cap(frm(y ~ x, family = asym_laplace(), data = d)))
}
set.seed(4000 + 7); d <- data.frame(x = rnorm(200), z = rnorm(200))
d$y <- 1 + 0.5 * d$x + exp(0.4 * d$z) * rnorm(200)
show("D7", cap(frm(bf(y ~ x, sigma ~ z), family = asym_laplace(), data = d)))
set.seed(5000 + 18); d <- data.frame(x = rnorm(200), z = rnorm(200))
d$y <- 1 + 0.5 * d$x + rnorm(200)
show("E18", cap(frm(bf(y ~ x, quantile ~ z), family = asym_laplace(), data = d)))
set.seed(5); d <- data.frame(x = rnorm(100)); d$y <- 1 + 0.5 * d$x + rnorm(100)
show("fixq-gauss5", cap(frm(bf(y ~ x, quantile = 0.5), family = asym_laplace(), data = d)))
for (s in 1:3) {
  set.seed(20 + s); n <- 60
  d <- data.frame(x = sample(0:2, n, TRUE)); d$y <- rpois(n, 3 + d$x)
  show(paste0("ties-fixq-s", s), cap(frm(bf(y ~ x, quantile = 0.5), family = asym_laplace(), data = d)))
  show(paste0("ties-estq-s", s), cap(frm(y ~ x, family = asym_laplace(), data = d)))
}
for (s in 1:6) {
  set.seed(40 + s); n <- 80
  d <- data.frame(x = sample(0:1, n, TRUE), g = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- rpois(n, 2 + d$x) * ifelse(d$g == "a", 1, 2)
  show(paste0("ties-sigmag-s", s), cap(frm(bf(y ~ x, sigma ~ g, quantile = 0.5), family = asym_laplace(), data = d)))
}
set.seed(31); d <- data.frame(y = round(rnorm(10), 3))
show("even-median", cap(frm(bf(y ~ 1, quantile = 0.5), family = asym_laplace(), data = d)))
set.seed(5); d <- data.frame(x = rnorm(100)); d$y <- 1 + 0.5 * d$x + rnorm(100)
f <- suppressWarnings(frm(y ~ x, family = asym_laplace(), data = d))
r <- cap(confint(f, parm = "x", method = "profile"))
cat("profile:", class(r$v)[1], substr(conditionMessage(r$v), 1, 120), "\n")
cat("DONE\n")
