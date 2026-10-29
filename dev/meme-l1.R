# Lane meme: L1 identity of frmtmb's me() posterior density, every class
# with a proper prior (meanme, sdme and corme included), against brms
# 2.23.0's Stan log_prob on the unconstrained scale, at several points,
# value and gradient. Output: dev/meme-out/l1.txt
#   Rscript dev/meme-l1.R > dev/meme-out/l1.txt
# The same check is gated test "row 23b" of test-brms-likelihood.R; this
# script prints the numbers.
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-meme/dev/stan-cache",
           NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
tt <- "C:/Users/adf44/source/r/frmtmb-wt-meme/tests/testthat"
te <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-me-lp.R")) {
  sys.source(file.path(tt, h), envir = te)
}
options(frmtmb.brms_lp_report = TRUE)
for (cs in te$me_lp_cases()) {
  cat("====", cs$label, "\n")
  r <- te$me_lp_points(cs, npts = 4, seed = 23)
  print(r$table, digits = 12)
  cat(sprintf("spread(const)/scale %.3g  max grad rel %.3g  max |Lnum - Lana| %.3g\n",
              r$spread, r$grad_rel, r$L_err))
}
