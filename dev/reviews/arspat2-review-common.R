# Review of lane arspat2: library setup. REVLIB="" gives the base build
# (rellib-r21) alone; default the lane library first.
lib <- Sys.getenv("REVLIB", "C:/Users/adf44/source/r/wt-arspat2-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
cat("lib:", dirname(dirname(system.file(package = "frmtmb"))),
    format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
`%||%` <- function(a, b) if (is.null(a)) b else a
rook <- function(k) {
  g <- expand.grid(r = seq_len(k), c = seq_len(k))
  n <- nrow(g)
  outer(seq_len(n), seq_len(n), function(i, j) {
    as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
  })
}
# objective of fit f against fit g at outer points near f's optimum, nu
# set through `inu` (log(nu - 1)) when given
cmp_pts <- function(f, g, inu = integer(0), nus = c(3, 6, 21),
                    map = identity) {
  out <- numeric(0)
  for (k in 0:2) {
    p <- f$opt$par + 0.05 * k * seq_along(f$opt$par) / length(f$opt$par)
    if (length(inu)) p[inu] <- log(nus[k + 1L] - 1)
    a <- -f$obj$fn(p)
    b <- -g$obj$fn(map(p))
    out <- c(out, a - b)
  }
  out
}
