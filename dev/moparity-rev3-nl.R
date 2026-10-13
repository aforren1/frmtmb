# Reviewer, final check (r3): nonlinear parameters named like specials,
# on data with no column of those names, lane against base.
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(3); n <- 300
d <- data.frame(z = rnorm(n)); d$ynl <- 2 * exp(-0.5 * abs(d$z)) + rnorm(n, 0, 0.2)
ll <- function(fo) {
  r <- tryCatch(suppressWarnings(suppressMessages(frm(fo, family = gaussian(), data = d))), error = function(e) e)
  if (inherits(r, "error")) paste("ERROR", substr(conditionMessage(r), 1, 110)) else format(as.numeric(logLik(r)), digits = 12)
}
cat("plain k:", ll(bf(ynl ~ a * exp(-k * abs(z)), a ~ 1, k ~ 1, nl = TRUE)), "\n")
for (p in c("s", "gp", "mo", "ar", "mi", "cs", "me")) {
  fo <- bf(as.formula(paste0("ynl ~ a * exp(-", p, " * abs(z))")), a ~ 1,
           as.formula(paste0(p, " ~ 1")), nl = TRUE)
  cat("nl parameter", p, ":", ll(fo), "\n")
}
