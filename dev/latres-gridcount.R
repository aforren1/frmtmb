# Who refuses: the cells brms builds and frm() refuses, by family, from the grid.
fs <- list.files("dev/latres-out/grid-after", "^grid-.*tsv$", full.names = TRUE)
g <- do.call(rbind, lapply(fs, function(f) {
  x <- readLines(f, warn = FALSE); m <- strsplit(x[-1L], "\t", fixed = TRUE)
  m <- lapply(m, function(v) c(v, rep("", 10L - length(v))))
  d <- as.data.frame(do.call(rbind, m), stringsAsFactors = FALSE)
  names(d) <- strsplit(x[1L], "\t", fixed = TRUE)[[1L]]; d }))
print(names(g))
g$bb <- g$brms != "refuses"; g$fr <- g$frm == "refuses"
cat("both build", sum(g$bb & !g$fr), " both refuse", sum(!g$bb & g$fr), " brms builds frm refuses", sum(g$bb & g$fr), " brms refuses frm fits", sum(!g$bb & !g$fr), "\n")
print(table(g$family[g$bb & g$fr]))
print(table(g$term[g$bb & g$fr]))
print(g[!g$bb & !g$fr, c("family", "term")])
