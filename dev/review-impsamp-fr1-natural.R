# Fix round 1, items 3 and 4: the lane's natural search, my run, on 7 designs at
# data seeds 1:12 and draws 2, 4, 8, 16, 64. Full messages, every capped fit's
# moves and warnings, and any fit whose move is not small yet reports nothing.
# Each refusal that names a count is refit at that count; the correlated slope
# fits are read against my own exact marginal at the end.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
args <- as.integer(commandArgs(TRUE))
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/impsamp-natural.R")
i <- grep("^mk <- list\\(", src)
j <- grep("^n <- c\\(fires", src)
mk <- eval(parse(text = sub("^mk <- ", "", paste(src[i:(j - 1L)], collapse = "\n"))))
out <- sprintf("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out/fr1/natural-%d.rds", args[1])
rows <- list()
for (nm in names(mk)[args[1]]) for (s in 1:12) for (N in c(2L, 4L, 8L, 16L, 64L)) {
  m <- mk[[nm]](s)
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(frm(m$f, family = m$fam, data = m$d, importance = N), error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  e <- list(design = nm, s = s, N = N, warn = w)
  if (inherits(r, "error")) {
    e$msg <- conditionMessage(r)
  } else {
    im <- r$importance
    e[c("rounds", "moves", "grads", "capped", "ess", "mcse", "ll", "conv", "par")] <-
      list(im$rounds, im$moves, im$grads, im$capped, im$ess_min, im$mcse,
           as.numeric(logLik(r)), r$opt$convergence, unname(r$opt$par))
  }
  rows[[length(rows) + 1L]] <- e
}
saveRDS(rows, out)
cat("done", args[1], length(rows), "\n")
