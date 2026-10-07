# Lane optima, item 4c: is a multi-start option worth having for an
# ordinal mixture? Per seed (two-class latent data of
# dev/optima-csmix.R, n = 300, no cs()): the default fit, then 10
# starts jittered about the default's start (every outer coordinate by
# N(0, 0.5), the generator frm_allfit()-style restarts would use). For
# the default and for the best of the 10: logLik, optimizer code, and
# whether the degenerate-component warning fired.
#   Rscript dev/optima-ordmix-starts.R lane|base seeds out.tsv [family]
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
seeds <- eval(parse(text = args[2]))
out <- args[3]
famname <- if (length(args) > 3) args[4] else "cum_cum"
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
fam <- function() {
  switch(famname,
    cum_cum = mixture(cumulative(), cumulative()),
    cum_sratio = mixture(cumulative(), sratio()))
}
fit1 <- function(d, start = NULL) {
  w <- character()
  f <- tryCatch(withCallingHandlers(
    frm(bf(y ~ x + z), family = fam(), data = d, start = start),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }), error = function(e) NULL)
  if (is.null(f)) return(list(ll = NA, code = NA, degen = NA))
  list(ll = as.numeric(logLik(f)), code = f$opt$convergence,
       degen = any(grepl("^mixture[(].*component", w)), fit = f)
}
rows <- list()
for (s in seeds) {
  d <- mk(s, mix = TRUE)
  f0 <- fit1(d)
  tpl <- frm(bf(y ~ x + z), family = fam(), data = d,
             dry_run = "objective")$estimates
  set.seed(1000 + s)
  best <- list(ll = -Inf)
  for (i in 1:10) {
    st <- lapply(tpl, function(v) v + stats::rnorm(length(v), 0, 0.5))
    # disc is held at 1 (mapped), so its slot keeps its value
    dk <- grepl("^disc", names(tpl$betad))
    st$betad[dk] <- tpl$betad[dk]
    r <- fit1(d, st)
    if (isTRUE(r$ll > best$ll)) best <- r
  }
  rows[[length(rows) + 1L]] <- data.frame(
    seed = s, ll0 = f0$ll, code0 = f0$code, degen0 = f0$degen,
    ll_best = best$ll, code_best = best$code, degen_best = best$degen)
}
X <- do.call(rbind, rows)
X$gain <- X$ll_best - X$ll0
utils::write.table(X, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(X)
g <- X$gain > 1e-3
cat("seeds", nrow(X), "; best of 10 above the default by > 1e-3:", sum(g),
    "; of those the best is degenerate:", sum(g & X$degen_best),
    "; default degenerate:", sum(X$degen0, na.rm = TRUE),
    "; best non-degenerate and above the default:",
    sum(g & !X$degen_best), "\n")
