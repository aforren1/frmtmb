fs <- list.files("dev/latres-out/r1", "tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.delim(f, stringsAsFactors = FALSE)))
for (fm in c("bb","beta","zoib","nb","gamma")) {
  n <- x[x$fam == fm & x$arm == "none", c("seed","ll","sderr","ar","ldisp")]
  m <- x[x$fam == fm & x$arm == "box" & x$cap == 15, c("seed","ll","sderr","ar","ldisp","code")]
  j <- merge(n, m, by = "seed", suffixes = c(".n", ".b"))
  j <- j[j$ldisp.n >= 14 & j$ldisp.n < 30, ]
  cat(fm, "seeds at cap 15 with an uncapped fit that stayed sane", nrow(j), " largest |ll diff|", max(abs(j$ll.n - j$ll.b)), " largest |sderr ratio-1|", max(abs(j$sderr.b / j$sderr.n - 1)), " |ar diff|", max(abs(j$ar.n-j$ar.b)), "\n")
}
