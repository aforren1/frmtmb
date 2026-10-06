# Reviewer of lane fixes, claim 4c: vcov_cluster(full = TRUE) on an
# ordinal fit, base against lane.
#   Rscript dev/fixes-rev-vcl.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
set.seed(31)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), cl = factor(rep(1:30, 10)))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
for (f in list(cumulative(), sratio())) {
  fit <- frm(y ~ x, family = f, data = d)
  r <- tryCatch(vcov_cluster(fit, cluster = d$cl, full = TRUE),
                error = function(e) conditionMessage(e))
  cat(f$family, "full = TRUE:", if (is.character(r)) r else rownames(r), "\n")
  r <- tryCatch(vcov_cluster(fit, cluster = d$cl),
                error = function(e) conditionMessage(e))
  cat(f$family, "full = FALSE:", if (is.character(r)) r else rownames(r), "\n")
}
fit <- frm(bf(y ~ x, disc ~ 0 + z), family = cumulative(), data = d)
r <- tryCatch(rownames(vcov_cluster(fit, cluster = d$cl, full = TRUE)),
              error = function(e) conditionMessage(e))
cat("disc model full = TRUE:", r, "\n")
