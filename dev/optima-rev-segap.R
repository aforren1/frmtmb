# Reviewer of lane optima, claim 8b: a mo() simplex with an unobserved
# interior category (four categories, the third never observed) has an
# exactly flat direction that is not along a coordinate. When the SE
# check is silent, what do the standard errors say?
#   Rscript dev/optima-rev-segap.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
se_mo_gap <- function(seed) {
  set.seed(seed)
  inc <- sample(c(0L, 1L, 3L), 100, TRUE)
  ls <- c(30, 60, 70, 75)[inc + 1L] + rnorm(100, sd = 7)
  data.frame(inc, ls, age = rnorm(100, mean = 40, sd = 10))
}
silent <- 0L
tot <- 0L
for (fo in list(ls ~ mo(inc), ls ~ mo(inc) + age, ls ~ mo(inc) * age)) {
  for (s in 1:6) {
    w <- character()
    f <- withCallingHandlers(frm(fo, data = se_mo_gap(s)),
                             warning = function(x) {
                               w <<- c(w, conditionMessage(x))
                               invokeRestart("muffleWarning")
                             })
    tot <- tot + 1L
    H <- f$obj$he(f$opt$par)
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    V <- tryCatch(vcov(f, full = TRUE), error = function(e) NULL)
    z <- grep("^zeta", rownames(V))
    sez <- sqrt(diag(V)[z])
    seb <- sqrt(diag(V)[-z])
    if (!length(w)) silent <- silent + 1L
    cat(sprintf(paste0("%-18s seed %d warn %d | H eig min/max %.3g | ",
                       "SE zeta %s | max SE other %.3g\n"),
                deparse(fo), s, length(w), min(abs(ev)) / max(abs(ev)),
                paste(format(sez, digits = 3), collapse = ","), max(seb)))
  }
}
cat("silent on", silent, "of", tot, "\n")
