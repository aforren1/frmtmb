# Lane optima, item 1: the plain monotonic fit, `ls ~ mo(income)` on
# brms_monotonic's data code, seeds 1 to 200, against its exact maximum
# (the best least-squares fit over both signs and every support of the
# three increments, which is exact by the KKT conditions), with the
# objective evaluations each fit spent.
#   Rscript dev/optima-mo-main.R base|lane seeds out.tsv
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
seeds <- eval(parse(text = args[2]))
out <- args[3]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-mo-study.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^# exact maximum", src) - 1L)]))
exact_main <- function(d) {
  code <- as.integer(d$income) - 1L
  S <- sapply(1:3, function(k) as.numeric(code >= k))
  n <- nrow(d)
  best <- -Inf
  for (s in c(-1, 1)) for (m in 0:7) {
    A <- which(bitwAnd(m, c(1L, 2L, 4L)) > 0)
    X <- cbind(1, S[, A, drop = FALSE])
    fit <- stats::lm.fit(X, d$ls)
    cf <- fit$coefficients[-1]
    if (length(A) && any(s * cf < 0)) next
    ll <- -(n / 2) * (log(2 * pi * sum(fit$residuals^2) / n) + 1)
    best <- max(best, ll)
  }
  best
}
rows <- list()
for (s in seeds) {
  d <- mk(s)
  f <- suppressWarnings(frm(ls ~ mo(income), data = d))
  rows[[length(rows) + 1L]] <- data.frame(
    seed = s, ll_fit = as.numeric(logLik(f)), ll_exact = exact_main(d),
    code = f$opt$convergence, evals = f$opt$evals,
    searched = !is.null(f$opt$mo_search))
}
X <- do.call(rbind, rows)
X$gap <- X$ll_exact - X$ll_fit
utils::write.table(X, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("seeds", nrow(X), "; gap <= 1e-6:", sum(X$gap <= 1e-6), "; max gap",
    format(max(X$gap), digits = 4), "; min gap", format(min(X$gap), digits = 4),
    "; codes", paste(names(table(X$code)), table(X$code), sep = "=",
                     collapse = " "), "; searched", sum(X$searched),
    "; evaluations total", sum(X$evals), "median", stats::median(X$evals),
    "\n")
