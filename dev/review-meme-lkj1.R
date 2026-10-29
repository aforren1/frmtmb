# Review of lane meme: does lkj(1) change frm()'s MAP for cor and rescor
# as it now does for corme? Rscript dev/review-meme-lkj1.R
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-meme/dev/review-meme-shapes.R")
d <- rv_data()
ll <- function(f, pr = NULL) {
  fit <- suppressWarnings(suppressMessages(frm(f, data = d, prior = pr)))
  sprintf("%.10f", as.numeric(logLik(fit)))
}
f1 <- bf(y ~ x + (1 + x | g))
cat("cor:    none", ll(f1), " lkj(1)", ll(f1, set_prior("lkj(1)", class = "cor")),
    "\n")
f2 <- mvbf(bf(y ~ x), bf(y2 ~ x)) + set_rescor(TRUE)
cat("rescor: none", ll(f2), " lkj(1)",
    ll(f2, set_prior("lkj(1)", class = "rescor")), "\n")
f3 <- bf(y ~ me(x, sx) + me(z, sz))
cat("corme:  none", ll(f3), " lkj(1)",
    ll(f3, set_prior("lkj(1)", class = "corme")), "\n")
cat("brms table rows for cor (as_priorlist):\n")
bt <- suppressMessages(brms::default_prior(brms::bf(y ~ x + (1 + x | g)),
                                           data = d))
cat("cor:    brms table", ll(f1, bt[bt$class != "sigma" &
                                       bt$class != "Intercept" &
                                       bt$class != "sd", ]), "\n")
