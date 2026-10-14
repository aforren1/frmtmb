# Summary of dev/relrev071-worst.R logs, per arm and case.
#   Rscript dev/relrev071-worst-sum.R <dir>
d <- commandArgs(trailingOnly = TRUE)[1]
for (f in list.files(d, "\\.txt$", full.names = TRUE)) {
  L <- readLines(f)[-1]
  L <- L[grepl("^deg", L)]
  fld <- function(k) trimws(vapply(strsplit(L, " \\| "), `[`, "", k))
  case <- sub(" +s *[0-9]+$", "", fld(1))
  num <- function(k, pre) as.numeric(sub(pre, "", fld(k)))
  x <- data.frame(case = case, ls = fld(3), code = num(4, "code "),
                  dll = num(6, "dLL "), rel = num(7, "maxrelSE "),
                  route = fld(8), lost = fld(9), warn = fld(10))
  cat("==", basename(f), "\n")
  for (cs in unique(x$case)) {
    y <- x[x$case == cs, ]
    cat(sprintf(paste0("  %-16s n %2d ls %2d code!=0 %d warned %d",
                       " min dLL %+.2e max relSE %.2e route %2d",
                       " lost boundary %2d\n"),
                cs, nrow(y), sum(y$ls == "ls TRUE"), sum(y$code != 0),
                sum(y$warn != "warn -"), min(y$dll), max(y$rel),
                sum(y$route == "route TRUE"),
                sum(grepl("boundary", y$lost))))
  }
}
