# Lane parity2, item 4337: where default_prior() fails on rr().
.libPaths(c("C:/Users/adf44/source/r/wt-parity2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(2)
drr <- data.frame(g = factor(rep(1:10, each = 6)),
                  tf = factor(rep(c("a", "b", "c"), 20)), x = rnorm(60))
drr$y <- rnorm(60)
r <- tryCatch(default_prior(bf(y ~ x + rr(0 + tf | g, d = 1)),
                            family = gaussian(), data = drr),
              error = function(e) {
                print(sys.calls())
                e
              })
print(r)
options(error = NULL)
withCallingHandlers(
  default_prior(bf(y ~ x + rr(0 + tf | g, d = 1)), family = gaussian(),
                data = drr),
  error = function(e) print(rev(sys.calls())[1:12]))
