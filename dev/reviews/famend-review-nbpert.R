# Review of lane famend: which perturbation of the nb/nb2 assertion of
# test-family-fitend.R fails, and why, under whichever R runs it.
#   Rscript dev/reviews/famend-review-nbpert.R
.libPaths(c("C:/Users/adf44/source/r/wt-famend-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "| R", R.home(), "\n")
conds <- function(expr) {
  w <- list(); m <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w[[length(w) + 1L]] <<- x; invokeRestart("muffleWarning")
  }, message = function(x) {
    m[[length(m) + 1L]] <<- x; invokeRestart("muffleMessage")
  })
  list(value = v, warnings = w, messages = m)
}
has <- function(l, cls) any(vapply(l, function(x) inherits(x, cls), NA))
set.seed(4)
dn <- data.frame(x = rnorm(300))
dn$y <- rpois(300, exp(1 + 0.3 * dn$x))
for (k in 0:20) {
  d <- dn
  d$x <- d$x * (1 + 2e-15 * k)
  r1 <- conds(frm(y ~ x, family = negbinomial(), data = d))
  r2 <- conds(frm(y ~ x, family = brmsfamily("negbinomial2"), data = d))
  lim <- function(r) has(r$messages, "frmtmb_limit_fit") ||
    has(r$messages, "frmtmb_boundary_fit")
  l1 <- as.numeric(logLik(r1$value))
  l2 <- as.numeric(logLik(r2$value))
  cat(sprintf(paste0("k %2d nb: code %d lim %s ll %.10f last %.4f | nb2: code ",
                     "%d lim %s ll %.10f last %.4f | rel %.3g%s\n"), k,
              r1$value$opt$convergence, lim(r1), l1,
              tail(r1$value$opt$par, 1), r2$value$opt$convergence, lim(r2), l2,
              tail(r2$value$opt$par, 1), abs(l1 - l2) / abs(l1),
              if (!(lim(r1) && lim(r2) && abs(l1 - l2) / abs(l1) < 1e-8))
                "  <== FAILS" else ""))
  if (!(lim(r1) && lim(r2))) {
    for (r in list(r1, r2)) {
      for (x in r$warnings) cat("   W:", substr(conditionMessage(x), 1, 200), "\n")
      for (x in r$messages) cat("   M:", substr(conditionMessage(x), 1, 200), "\n")
    }
  }
}
cat("END\n")
