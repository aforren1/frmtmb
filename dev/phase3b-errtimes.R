# Punch round 1, minor 3: every first-pass eam recovery error, its
# kind and the time its RDS was written. Output: dev/phase3b-log/errtimes.txt
f <- list.files("dev/phase3b-log/recov", "^(cens|cont|contdef|contfix|collapse)-",
                full.names = TRUE)
rows <- lapply(f, function(p) {
  x <- readRDS(p)
  if (is.null(x$error)) return(NULL)
  data.frame(file = basename(p),
             kind = if (grepl("bad_alloc", x$error)) "bad_alloc" else
               if (grepl("NA/NaN gradient", x$error)) "NaN gradient" else x$error,
             written = format(file.mtime(p), "%m-%d %H:%M"))
})
tab <- do.call(rbind, rows)
tab <- tab[order(tab$written), ]
out <- c(capture.output(print(tab, row.names = FALSE)), "",
         capture.output(print(table(tab$kind))))
writeLines(out, "dev/phase3b-log/errtimes.txt")
cat(out, sep = "\n")
