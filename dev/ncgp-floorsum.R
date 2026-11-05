# The reviewer's table of M1 on the warm-started build, from dev/review-ncgp-floor.R's logs.
x <- unlist(lapply(c("dev/ncgp-out/m1/floor-40.txt", "dev/ncgp-out/m1/floor-100.txt"), readLines, warn = FALSE))
x <- grep("^RF ", x, value = TRUE)
num <- function(p) suppressWarnings(as.numeric(sub(paste0(".*", p, " +([^ ]+).*"), "\\1", x)))
sg <- num("frm sigma/ysc"); sd <- num("sd/ysc"); conv <- num(" conv"); nw <- num("nwarn")
gap <- suppressWarnings(as.numeric(sub(".*gap +([^ ]+).*", "\\1", x)))
cat(sprintf("rows %d, a sample line: %s\n", length(x), substr(x[1], 1, 170)))
r <- sg / sd
cat(sprintf("FS worst gap %.3g; short>1e-3 %d; short and silent %d\n", min(gap, na.rm = TRUE), sum(gap < -1e-3, na.rm = TRUE),
            sum(gap < -1e-3 & nw == 0, na.rm = TRUE)))
for (b in list(c(0, 1e-6), c(1e-6, 1e-5), c(1e-5, 1e-4), c(1e-4, Inf))) {
  k <- which(r >= b[1] & r < b[2])
  cat(sprintf("FS sigma/sdgp %.0e to %.0e rows %3d warned %3d conv1 %3d short>1e-3 %3d short silent %d worst %.3g\n", b[1], b[2], length(k),
              sum(nw[k] > 0), sum(conv[k] != 0), sum(gap[k] < -1e-3, na.rm = TRUE), sum(gap[k] < -1e-3 & nw[k] == 0, na.rm = TRUE), min(gap[k], na.rm = TRUE)))
}
