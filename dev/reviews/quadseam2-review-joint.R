# Review of lane quadseam2, attack 7: the merge with selimit2. Its joint
# ends (se_limit_joint_dirs()) read `fn`, which this lane points at rules
# taped at each point under one-block quadrature; the combination was run
# by neither lane. A treatment-coded poisson whose reference level is all
# zero (intercept and the other levels run off together), and a
# bernoulli `zi`-free analogue, by Laplace and by quadrature, seeds 1-6.
#   Rscript quadseam2-review-joint.R <libpath>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
verdict <- function(fo, fam, d, quad) {
  w <- character()
  f <- withCallingHandlers(
    frm(fo, family = fam, data = d, quadrature = quad),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    }, message = function(x) {
      w <<- c(w, paste0("m:", class(x)[1])); invokeRestart("muffleMessage")
    })
  se <- withCallingHandlers(sqrt(diag(vcov(f))),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    }, message = function(x) {
      w <<- c(w, paste0("m:", class(x)[1])); invokeRestart("muffleMessage")
    })
  sprintf("nonfinite[%s] w[%s] se[%s]",
          paste(names(se)[!is.finite(se)], collapse = ","),
          paste(sort(unique(w)), collapse = ","),
          paste(format(se, digits = 3), collapse = ","))
}
for (s in 1:6) {
  set.seed(s)
  d <- data.frame(h = factor(rep(1:20, each = 9)),
                  f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
  d$y <- rpois(180, exp(0.5 + 0.3 * d$x + rnorm(20, 0, 0.5)[d$h]))
  d$y[d$f == "a"] <- 0L
  for (q in c(FALSE, TRUE)) {
    cat(sprintf("POISREF %d %s %s\n", s, if (q) "QUAD" else "LAP",
                verdict(y ~ f + x + (1 | h), poisson(), d, q)))
  }
}
cat("END\n")
