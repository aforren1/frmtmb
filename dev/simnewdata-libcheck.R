# After the machine crash of 2026-09-24 09:24: does the lane library
# still hold what the worktree source says, and do the three libraries
# have hollow packages? Compares every function the lane's changed R
# files define, source against installed namespace, by deparsed code.
#   Rscript dev/simnewdata-libcheck.R
source("dev/simnewdata-prelude.R")
libs <- c("C:/Users/adf44/source/r/simnewdata-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
for (L in libs) {
  d <- list.dirs(L, recursive = FALSE, full.names = FALSE)
  h <- d[!file.exists(file.path(L, d, "DESCRIPTION"))]
  cat(L, length(d), "dirs,", length(h), "hollow", head(h), "\n")
}
check <- function(pkg, files) {
  ns <- asNamespace(pkg)
  n_ok <- 0L
  bad <- character(0)
  for (f in files) {
    e <- new.env(parent = ns)
    exprs <- parse(f, keep.source = FALSE)
    for (ex in exprs) {
      if (is.call(ex) && identical(ex[[1L]], as.name("<-")) &&
          is.name(ex[[2L]]) && is.call(ex[[3L]]) &&
          identical(ex[[3L]][[1L]], as.name("function"))) {
        nm <- as.character(ex[[2L]])
        src <- eval(ex[[3L]], e)
        got <- get0(nm, envir = ns, inherits = FALSE)
        same <- !is.null(got) &&
          identical(deparse(src), deparse(got))
        if (same) n_ok <- n_ok + 1L else bad <- c(bad, paste(f, nm))
      }
    }
  }
  cat(pkg, "from", find.package(pkg), ":", n_ok, "functions match;",
      length(bad), "differ", if (length(bad)) paste(bad, collapse = "; "),
      "\n")
}
suppressMessages(library(frmtmb))
check("frmtmb", c("R/predict.R", "R/simulate-newdata.R", "R/re-formula.R",
                  "R/conditional-effects.R", "R/bootstrap.R"))
check("frmtmb.sample", "extensions/frmtmb.sample/R/methods-draws.R")
check("frmtmb.learn", c("extensions/frmtmb.learn/R/family.R"))
# and a fit that produces a number
set.seed(1)
d <- data.frame(x = rnorm(200), g = factor(rep(1:10, 20)))
d$y <- rnorm(200, 1 + d$x + rnorm(10)[d$g])
f <- frm(bf(y ~ x + (1 | g)), data = d)
cat("gaussian (1 | g) logLik", format(as.numeric(logLik(f)), digits = 12),
    "\n")
