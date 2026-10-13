# rel070 seflags: run one of the dev/rel070-seflags-* scripts unchanged,
# except that cap() also prints every warning it catches from a frm()
# call ("FITW:" lines), so the fit-time warnings the script only counts
# can be read. Everything else the script prints is as it was.
#   Rscript rel070-seflags-capwarn.R <script.R> <script args...>
a <- commandArgs(trailingOnly = TRUE)
script <- a[1]
sargs <- a[-1]
commandArgs <- function(trailingOnly = FALSE) {
  if (trailingOnly) sargs else c(base::commandArgs(FALSE), sargs)
}
src <- readLines(script)
i <- grep("^source[(]", src)[1]
eval(parse(text = src[i]), envir = globalenv())
cap0 <- cap
cap <- function(expr) {
  e <- substitute(expr)
  r <- eval(substitute(cap0(E), list(E = e)), parent.frame())
  lab <- paste(deparse(e, width.cutoff = 500L), collapse = " ")
  if (grepl("^frm[(]", lab)) {
    for (w in r$w) cat("    FITW:", substr(gsub("\n", " ", w), 1, 260), "\n")
    if (!inherits(r$v, "rev_err") && !is.null(r$v$opt$polish)) {
      cat("    POLISH:", paste(names(r$v$opt$polish),
                                signif(r$v$opt$polish, 4), collapse = " "),
          "\n")
    }
    if (!inherits(r$v, "rev_err")) {
      u <- asNamespace("frmtmb")$predictor_units
      if (!is.null(u)) {
        uu <- u(r$v$frame)
        cat("    UNITS:", if (length(uu)) paste(names(uu), uu) else "-", "\n")
      }
    }
  }
  r
}
eval(parse(text = src[-seq_len(i)]), envir = globalenv())
