# The table of dev/ncgp-nuggetrate.R: by how far sigma reads low, how often the nugget
# warning fires, and how often any warning does.
x <- readLines("dev/ncgp-out/nuggetrate.txt", warn = FALSE)
x <- grep("^NR ", x, value = TRUE)
num <- function(k) as.numeric(sub(paste0(".* ", k, " ([^ ]+).*"), "\\1", x))
ratio <- suppressWarnings(num("sigma/s"))
sgsd <- num("sigma/sdgp")
s <- num("s")
nw <- num("nugget-warning")
w <- num("warns")
band <- cut(ratio, c(0, 0.2, 0.5, 0.8, 1.25, Inf), right = FALSE)
band <- addNA(band)
for (b in levels(band)) {
  k <- which(band == b | (is.na(b) & is.na(band)))
  cat(sprintf("NRT sigma / true %-12s rows %3d nugget warning %3d any warning %3d\n",
              if (is.na(b)) "noise-free" else b, length(k), sum(nw[k] > 0), sum(w[k] > 0)))
}
for (lo in list(c(0, 1e-7), c(1e-7, 1e-6), c(1e-6, 3e-6), c(3e-6, 1e-5), c(1e-5, 1e-4), c(1e-4, Inf))) {
  k <- which(sgsd >= lo[1] & sgsd < lo[2])
  cat(sprintf("NRU sigma/sdgp %.0e to %.0e rows %3d nugget warning %3d any warning %3d\n",
              lo[1], lo[2], length(k), sum(nw[k] > 0), sum(w[k] > 0)))
}
for (v in sort(unique(s))) {
  k <- which(s == v)
  cat(sprintf("NRV true s %.0e rows %3d median sigma/true %s nugget warning %3d any %3d\n", v, length(k),
              format(stats::median(ratio[k]), digits = 3), sum(nw[k] > 0), sum(w[k] > 0)))
}
