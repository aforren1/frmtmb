# Lane selimit2, backlog 3389 and 3575: what an asym_laplace() fit says
# on 0.82.0 and the lane. dev/cdfs-dbg14.R's data (seed 43, 160 rows
# drawn from the family, quantile 0.3, uncensored), test-se-limit.R's
# (seed 5, gaussian data), and a zero_inflated_asym_laplace fit. Per fit:
# the code and message, the gradient, the lost reasons, every warning and
# message at fit time, and what summary() and fitted() say.
#   Rscript dev/selimit2-asym.R <base|lane>
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, dirname(find.package("frmtmb")), "\n")
cap <- function(expr) {
  w <- character(); m <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(z) {w <<- c(w, conditionMessage(z))
      invokeRestart("muffleWarning")},
    message = function(z) {m <<- c(m, conditionMessage(z))
      invokeRestart("muffleMessage")})
  list(v = v, w = w, m = m)
}
say <- function(lab, r) {
  for (x in r$w) cat("  ", lab, " WARN: ", gsub("\n", " ", x), "\n", sep = "")
  for (x in r$m) cat("  ", lab, " MSG: ", gsub("\n", " ", x), "\n", sep = "")
}
cases <- list()
set.seed(43); n <- 160; x <- rnorm(n)
y <- brms::rskew_normal(n, 1 + 0.5 * x, 1, 4)
y <- brms::rexgaussian(n, 2 + 0.5 * x, 0.5, 1)
y <- brms::rasym_laplace(n, 1 + 0.5 * x, 1, 0.3)
cases$dbg14 <- list(fo = y ~ x, fam = asym_laplace(),
                    d = data.frame(y = y, x = x))
set.seed(5)
d <- data.frame(x = stats::rnorm(100))
d$y <- 1 + 0.5 * d$x + stats::rnorm(100)
cases$gauss5 <- list(fo = y ~ x, fam = asym_laplace(), d = d)
cases$fixq <- list(fo = bf(y ~ x, quantile = 0.5), fam = asym_laplace(),
                   d = d)
set.seed(9)
dz <- data.frame(x = stats::rnorm(200))
dz$y <- ifelse(stats::runif(200) < 0.2, 0,
               brms::rasym_laplace(200, 1 + 0.5 * dz$x, 1, 0.5))
cases$zi <- list(fo = y ~ x, fam = zero_inflated_asym_laplace(), d = dz)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  r <- cap(frm(cs$fo, family = cs$fam, data = cs$d))
  if (inherits(r$v, "error")) {
    cat(nm, "ERROR", conditionMessage(r$v), "\n")
    next
  }
  f <- r$v
  g <- tryCatch(max(abs(f$obj$gr(f$opt$par))), error = function(e) NA)
  lost <- ns$sdr_of(f)$se_lost
  cat(sprintf("%s | code %d (%s) | max |grad| %.2e | lost %s\n", nm,
              f$opt$convergence, f$opt$message, g,
              if (length(lost)) paste(names(lost), lost, sep = ":",
                                      collapse = ",") else "none"))
  say("fit", r)
  s <- cap(summary(f))
  say("summary", s)
  fp <- cap(stats::fitted(f))
  if (!inherits(fp$v, "error")) {
    cat(sprintf("  fitted: %d of %d SEs NaN\n",
                sum(is.nan(fp$v[, "Est.Error"])), nrow(fp$v)))
  }
  say("fitted", fp)
}
