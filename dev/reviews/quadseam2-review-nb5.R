# Review of lane quadseam2, attack 3: the one verdict the lane changed in
# quadseam2-review-touch.R, negbinomial true sd 0 seed 5 (a new
# frmtmb_limit_fit message under quadrature). Does the Laplace fit of
# the same data say the same? Messages, the shape and the SEs.
#   Rscript quadseam2-review-nb5.R <base|lib>
a <- commandArgs(TRUE)
if (a[1] == "lib") a[1] <- "C:/Users/adf44/source/r/wt-quadseam2-lib"
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(5)
d <- data.frame(g = factor(rep(1:15, each = 4)), x = rnorm(60))
u <- rnorm(15, 0, 0)[d$g]
d$y <- rnbinom(60, mu = exp(0.3 + 0.5 * d$x + u), size = 3)
for (q in c(FALSE, TRUE)) {
  m <- character()
  f <- withCallingHandlers(frm(y ~ x + (1 | g), family = negbinomial(),
                               data = d, quadrature = q),
    message = function(x) {
      m <<- c(m, paste(class(x)[1], substr(gsub("\\s+", " ",
                                                  conditionMessage(x)), 1, 160)))
      invokeRestart("muffleMessage")
    })
  se <- withCallingHandlers(sqrt(diag(vcov(f))), message = function(x) {
    m <<- c(m, paste(class(x)[1], substr(gsub("\\s+", " ",
                                                conditionMessage(x)), 1, 160)))
    invokeRestart("muffleMessage")
  })
  cat(sprintf("QUAD %s ll %.9f betad %s\n  se %s\n", q,
              as.numeric(logLik(f)),
              paste(format(f$estimates$betad, digits = 5), collapse = ","),
              paste(names(se), format(se, digits = 3), collapse = " ")))
  for (x in unique(m)) cat("  ", x, "\n")
}
cat("END\n")
