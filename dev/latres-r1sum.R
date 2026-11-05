# Summary of dev/latres-r1.R's arms.
fs <- list.files("dev/latres-out/r1", "tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.delim(f, stringsAsFactors = FALSE)))
x$over <- x$ll - x$ref
tab <- do.call(rbind, lapply(split(x, list(x$fam, x$arm, x$cap), drop = TRUE), function(z) {
  data.frame(fam = z$fam[1], arm = z$arm[1], cap = z$cap[1], fits = sum(is.finite(z$ll)),
             err = sum(nzchar(z$err) & !is.na(z$err)), code1 = sum(z$code != 0, na.rm = TRUE),
             ll_max = max(z$ll, na.rm = TRUE), ll_med = median(z$ll, na.rm = TRUE),
             ldisp_max = max(z$ldisp, na.rm = TRUE), over_med = median(z$over, na.rm = TRUE),
             over_max = suppressWarnings(max(abs(z$over), na.rm = TRUE)))
}))
tab <- tab[order(tab$fam, tab$arm, tab$cap), ]
print(format(tab, digits = 4), row.names = FALSE)
# each seed against the no-cap arm: how far a capped arm moved a fit that never
# reached the cap
for (fm in unique(x$fam)) {
  n <- x[x$fam == fm & x$arm == "none", c("seed", "ll", "sderr", "ldisp")]
  for (a in c("nan", "box")) for (k in c(15, 25)) {
    m <- x[x$fam == fm & x$arm == a & x$cap == k, c("seed", "ll", "sderr", "ldisp")]
    j <- merge(n, m, by = "seed", suffixes = c(".n", ".m"))
    j <- j[is.finite(j$ldisp.n) & j$ldisp.n < k - 1, ]
    if (nrow(j)) cat(sprintf("%-8s %s %d: seeds under the cap %d, largest |ll change| %.3g\n", fm, a, k, nrow(j),
                             max(abs(j$ll.n - j$ll.m), na.rm = TRUE)))
  }
}
