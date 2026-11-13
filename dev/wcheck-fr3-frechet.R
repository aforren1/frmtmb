# Lane wcheck, fix round 3: frechet's nu under its default link, logm1
# (nu = 1 + exp(eta), dnu/deta = nu - 1, so the MAP objective's
# Jacobian power is 1 as nu grows). A frechet density narrows like
# nu^-1 (rate 1), so at an excess E the likelihood grows like nu^E, and
# student_t(df) on nu (p = df + 1) gives the MAP objective nu^(E + 1 -
# p). Each arm: the check's verdict on both routes, then the slope of
# the profile along log nu (the rest re-optimized) against that power.
#   Rscript dev/wcheck-fr3-frechet.R <lib|base>
# Data: dev/reviews/wcheck-review-rates.R's frechet column (seed
# 20261111, n = 80, nu = 4), weights 1 + E / n.
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n))
mu <- exp(0.5 + 0.3 * d$x + rnorm(n, 0, 0.3))
d$yf <- mu * (-log(runif(n)))^(-1 / 4)
E <- 6
d$w <- 1 + E / n
fo <- bf(yf | weights(w) ~ x + (1 | obs))
fr <- q(frm(fo, family = frechet(), data = d, dry_run = "frame"))
cat("link of nu:", fr$linpreds[[frmtmb:::linpred_key(names(fr$spec$responses)[1], "nu")]]$link$name, "\n")
cw <- frmtmb:::check_weights_unbounded
for (df in c(5.5, 6.5)) {
  pr <- set_prior(sprintf("student_t(%g, 0, 2.5)", df), class = "nu")
  rt <- vapply(c("map", "sample"), function(route) {
    tryCatch({cw(fr$spec, fr, pr, exact = FALSE, route = route); "pass"},
             error = function(e) "REFUSED")
  }, "")
  fit <- tryCatch({
    f <- q(frm(fo, family = frechet(), data = d, prior = pr))
    sprintf("fit logLik %.6g conv %d", as.numeric(logLik(f)),
            f$opt$convergence)
  }, error = function(e) "REFUSED")
  ob <- q(frm(fo, family = frechet(), data = d, prior = pr,
              dry_run = "objective"))
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  k <- which(names(p) == "betad")[1L]
  free <- setdiff(seq_along(p), k)
  lx <- c(4, 6, 8, 10, 12)
  v <- vapply(log(exp(lx) - 1), function(t) {
    p[k] <- t
    fr2 <- function(r) {
      q2 <- p; q2[free] <- r
      v <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(v)) 1e10 else v
    }
    op <- tryCatch(nlminb(p[free], fr2), error = function(e) NULL)
    if (is.null(op)) return(NA_real_)
    if (op$objective < 1e9) p[free] <<- op$par
    -op$objective
  }, 0)
  cat(sprintf("frechet E %g, student_t(%g) on nu | map %s, sample %s | %s\n  slope %s | predicted %g\n",
              E, df, rt[["map"]], rt[["sample"]], fit,
              paste(format(diff(v) / diff(lx), digits = 3), collapse = " "),
              E + 1 - (df + 1)))
}
cat("DONE\n")
