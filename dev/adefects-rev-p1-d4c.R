source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages(library(frmtmb))
addr <- function(x) {
  strsplit(trimws(utils::capture.output(.Internal(inspect(x)))[1L]),
           "[[:space:]]+")[[1L]][1L]
}
set.seed(20260917)
d <- data.frame(g = factor(rep(1:6, each = 5)), x = rnorm(30))
d$y <- d$x + rnorm(30)
fit <- frm(y ~ x + (1 | g), d)
mf <- fit$data
# a copy that keeps the SAME terms environment: touch a column, which
# duplicates the list and its columns but not the attribute's env
cp <- mf
cp[[1L]][1L] <- cp[[1L]][1L]
cat("addresses differ            :", !identical(addr(mf), addr(cp)), "\n")
cat("identical(mf, cp)           :", identical(mf, cp), "\n")
cat("same terms env              :",
    identical(environment(attr(mf, "terms")),
              environment(attr(cp, "terms"))), "\n")
