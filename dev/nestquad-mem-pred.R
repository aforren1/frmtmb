# The final guard's prediction for each measured fit of
# dev/nestquad-mem3.log, from its logged tape bytes and stacked rows
# (arithmetic on the logged numbers, not a new measurement).
.libPaths(c("C:/Users/adf44/source/r/wt-nestquad-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
a <- commandArgs(TRUE)
x <- readLines(file.path("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev",
                        if (length(a)) a[1] else "nestquad-mem3.log"))
get <- function(l, k) {
  w <- strsplit(trimws(l), " +")[[1]]
  w[which(w == k)[1] + 1]
}
for (l in x) {
  key <- get(l, "key"); n <- as.numeric(get(l, "n"))
  delta <- as.numeric(get(l, "delta")); tape <- as.numeric(get(l, "tape"))
  plan <- list(row_inner = rep(1L, n))
  p <- ns$quad_nest_bytes(tape * 2^20, plan, c(11L, 15L), key) / 2^20
  pf <- ns$quad_nest_bytes(NULL, plan, c(11L, 15L), key) / 2^20
  cat(sprintf("PRED %-26s n %5d measured %6.0f MB predicted %6.0f MB ratio %.2f | without tape size %6.0f MB ratio %.2f\n",
              key, n, delta, p, p / delta, pf, pf / delta))
}
