# The whole warning a fit raises when nlminb stops at false convergence at true sigma
# 3e-4, on dev/ncgp-perm2.R's designs: what does the user read, and is it true?
#   Rscript dev/ncgp-fcwarn.R
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
for (n in c(40L, 100L, 200L)) for (seed in 1:5) {
  set.seed(seed)
  x <- round(stats::runif(n, 0, 6), 3)
  y <- sin(2 * x) + 3e-4 * stats::rnorm(n)
  w <- character()
  fit <- withCallingHandlers(
    frm(bf(y ~ gp(x)), data = data.frame(x = x, y = y), family = gaussian()),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  cat(sprintf("NCGPFC n %d seed %d conv %d warnings %d: %s\n", n, seed,
              fit$opt$convergence, length(w), paste(gsub("[\r\n]+", " ", w), collapse = " || ")))
}
