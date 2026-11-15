# Review of lane famend, item 3060: what reads the moved last.par.best
# after frechet_fit_check() on an edge fit. Lane build; "off" replaces the
# check by a no-op in this process.
#   Rscript dev/reviews/famend-review-frechet2.R
lanelib <- "C:/Users/adf44/source/r/wt-famend-lib"
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(lanelib, base, ulib))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
real_check <- get("frechet_fit_check", ns)
q <- function(expr) suppressMessages(suppressWarnings(expr))
set.seed(9031)
d <- data.frame(x = rnorm(300))
d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / 1.1) * rexp(300)^(-1 / 1.1)
out <- list()
for (on in c(TRUE, FALSE)) {
  assignInNamespace("frechet_fit_check",
                    if (on) real_check else function(fit, resp) invisible(NULL),
                    ns)
  f <- q(frm(y ~ x, family = frechet(), data = d))
  k <- if (on) "on" else "off"
  ci <- tryCatch(q(confint(f, parm = "x", method = "profile")),
                 error = function(e) conditionMessage(e))
  set.seed(1)
  pe <- tryCatch(q(predict(f, newdata = d[1:3, , drop = FALSE],
                           ndraws = 200)),
                 error = function(e) conditionMessage(e))
  set.seed(1)
  fl <- tryCatch(q(frm_linpred(f, newdata = d[1:3, , drop = FALSE],
                               se.fit = TRUE)),
                 error = function(e) conditionMessage(e))
  dg <- tryCatch(q(capture.output(diagnose(f))), error = function(e)
    conditionMessage(e))
  out[[k]] <- list(ci = ci, pe = pe, fl = fl, dg = dg,
                   sm = q(capture.output(summary(f))))
  cat("\n== check", k, "\n")
  print(ci)
  print(pe)
}
for (nm in names(out$on)) {
  cat(nm, "identical on/off:", identical(out$on[[nm]], out$off[[nm]]), "\n")
}
cat("\ndiagnose on:\n"); cat(out$on$dg, sep = "\n")
cat("\ndiagnose off:\n"); cat(out$off$dg, sep = "\n")
cat("\ndone\n")
