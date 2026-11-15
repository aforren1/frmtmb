# Review of lane famend, item 3075: does the profile interval the message
# points to work on these fits, and do refits repeat the message?
#   Rscript dev/reviews/famend-review-gevmsg.R <base|lane>
a <- commandArgs(trailingOnly = TRUE)
lanelib <- "C:/Users/adf44/source/r/wt-famend-lib"
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] == "lane") lanelib, base, ulib))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
cnt <- function(expr) {
  m <- character(); w <- character()
  v <- tryCatch(withCallingHandlers(expr,
    message = function(x) {
      m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
    }, warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }), error = function(e) structure(conditionMessage(e), class = "err"))
  cat("   messages", length(m), "( nonregular",
      sum(grepl("below -0.5", m)), ") warnings", length(w), "\n")
  for (x in unique(w)) cat("   W:", substr(gsub("\\s+", " ", x), 1, 160), "\n")
  v
}
rg <- function(n, mu, s, xi) mu + s * ((-log(runif(n)))^(-xi) - 1) / xi
for (xi0 in c(-0.63, -0.89)) {
  set.seed(4)
  d <- data.frame(y = rg(300, 0, 1, xi0), x = rnorm(300))
  cat("\n== xi0", xi0, "fit\n")
  f <- cnt(frm(y ~ 1, family = gen_extreme_value(), data = d))
  cat("== wald\n")
  print(cnt(confint(f)))
  cat("== profile\n")
  print(cnt(confint(f, method = "profile")))
  cat("== update(f, . ~ . + x)\n")
  invisible(cnt(update(f, . ~ . + x)))
  cat("== predicted xi: xi ~ x\n")
  invisible(cnt(frm(bf(y ~ 1, xi ~ x), family = gen_extreme_value(), data = d)))
}
cat("\ndone\n")
