# Short check of lane quadseam's fix round, item 3: the rebuilt runaway
# test's case (gaussian mi() with 3 rows missing, weights 1.5, (1 | obs),
# seed 20261111) on a given library: refused, or fitted, and how.
#   Rscript quadseam-short-wcheck.R <libpath>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$y[c(5, 17, 33)] <- NA
d$w <- 1.5
for (quad in c(TRUE, FALSE)) {
  w <- character()
  r <- tryCatch(withCallingHandlers(
    frm(bf(y | weights(w) + mi() ~ x + (1 | obs)), family = gaussian(),
        data = d, quadrature = quad),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")),
    error = function(e) e)
  if (inherits(r, "error")) {
    cat("QUAD", quad, "REFUSED:", substr(conditionMessage(r), 1, 200), "\n")
  } else {
    cat(sprintf("QUAD %s FIT logLik %.4f log sigma %.3f conv %d warn %s\n",
                quad, as.numeric(logLik(r)), r$estimates$betad[1],
                r$opt$convergence, paste(w, collapse = ",")))
  }
}
