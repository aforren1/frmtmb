# Reviewer of lane setier: an active upper bound on the top coefficient of
# a raw polynomial, degree 3 (well conditioned) and degree 5 (an
# eigenvalue 1e-11 of the largest): is the bound's verdict the same?
#   Rscript dev/setier-rev-bound.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, "\n")
for (deg in c(2, 3, 5)) {
  set.seed(5)
  d <- data.frame(x = runif(200, 1, 2))
  d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  top <- paste0("x", deg)
  fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + ")))
  f0 <- frm(fo, data = d)
  b <- fixef(f0)[top, "Estimate"]; s <- fixef(f0)[top, "Est.Error"]
  w <- character()
  f <- withCallingHandlers(frm(fo, data = d,
         prior = set_prior("", class = "b", coef = top, ub = b - s)),
       warning = function(x) {w <<- c(w, conditionMessage(x))
         invokeRestart("muffleWarning")})
  fe <- fixef(f)
  lost <- ns$sdr_of(f)$se_lost
  cat(sprintf("degree %d: %s at its bound %.4g (est %.4g) | SE %s | lost %s | %s\n",
              deg, top, b - s, fe[top, "Estimate"],
              paste(signif(fe[, "Est.Error"], 4), collapse = " "),
              if (length(lost)) paste(names(lost), lost) else "none",
              if (length(w)) substr(w[1], 1, 80) else "silent"))
}
