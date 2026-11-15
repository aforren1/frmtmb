# Short check of lane arspat2's fix round: library setup. SHORT_LIB=""
# gives the base build (rellib-r21) alone; default the lane library first.
lib <- Sys.getenv("SHORT_LIB", "C:/Users/adf44/source/r/rellib-r22")
.libPaths(c(if (nzchar(lib)) lib else "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
cat("lib:", dirname(system.file(package = "frmtmb")),
    format(packageVersion("frmtmb")), "| spatial:",
    dirname(system.file(package = "frmtmb.spatial")),
    "| margin slot:", !is.null(formals(frmtmb_register_autocor)$margin) ||
      "margin" %in% names(formals(frmtmb_register_autocor)), "\n")
ns <- asNamespace("frmtmb")
arm <- if (nzchar(lib)) "lane" else "base"
outdir <- "dev/rel083-log/kt/arspat2"
dir.create(outdir, showWarnings = FALSE)
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
rook <- function(k) {
  g <- expand.grid(r = seq_len(k), c = seq_len(k))
  n <- nrow(g)
  outer(seq_len(n), seq_len(n), function(i, j) {
    as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
  })
}
# outer points: the optimum moved by 0.05 k (scaled), k = 0, 1, 2
pts <- function(par, k) par + 0.05 * k * seq_along(par) / length(par)
reld <- function(a, b) abs(a - b) / abs(b)
peak_gb <- function() ps::ps_memory_info()[["peak_wset"]] / 2^30
