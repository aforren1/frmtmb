# The hmm_starts() prior test under three arms, one per process:
#   base   rellib-r18 core and frmtmb.latent
#   mixed  the lane's core with rellib-r18's frmtmb.latent (logLik()
#          is now the data part, and the old hmm_starts() reads it)
#   lane   the lane's core and frmtmb.latent
# Usage: Rscript dev/loglik-r80-latent-arms.R <arm>
arm <- commandArgs(TRUE)[1]
mine <- "C:/Users/adf44/source/r/wt-loglik-lib"
r18 <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(r18, user))
core_lib <- if (arm == "base") r18 else mine
lat_lib <- if (arm == "lane") mine else r18
suppressMessages({
  library(testthat)
  library(frmtmb, lib.loc = core_lib)
  library(frmtmb.latent, lib.loc = lat_lib)
})
cat("arm:", arm, " frmtmb", dirname(find.package("frmtmb")),
    " frmtmb.latent", dirname(find.package("frmtmb.latent")), "\n")
f <- "extensions/frmtmb.latent/tests/testthat/test-hmm-starts.R"
# only the prior test: the file's helpers, then that block
src <- parse(f, keep.source = FALSE)
env <- testthat::test_env("frmtmb.latent")
for (e in src) {
  is_test <- is.call(e) && identical(e[[1L]], as.name("test_that"))
  if (!is_test) {
    eval(e, env)
  } else if (grepl("under a prior the refits", e[[2L]], fixed = TRUE)) {
    res <- tryCatch(eval(e, env), error = function(err) err)
  }
}
cat("done\n")
