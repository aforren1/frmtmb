## Per-package totals of a dev/mocode-par.sh tier log, and every file
## with a fail, an error, a warning or no RESULT line.
##   Rscript dev/mocode-summary.R dev/mointerop-out/<tier>.log
f <- commandArgs(trailingOnly = TRUE)[1L]
x <- readLines(f)
cat("source:", f, "\n")
cat(grep(" ran ", x, value = TRUE), "\n")
r <- grep(" RESULT ", x, value = TRUE)
num <- function(k) as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", r))
d <- data.frame(pkg = sub(" .*", "", r),
                file = sub(".* RESULT ([^ ]+) .*", "\\1", r),
                pass = num("pass"), fail = num("fail"), error = num("error"),
                skip = num("skip"), warn = num("warn"))
a <- aggregate(cbind(files = 1L, pass, fail, error, skip, warn) ~ pkg,
               data = cbind(d, files = 1L), FUN = sum)
print(a, row.names = FALSE)
bad <- d[d$fail > 0 | d$error > 0 | d$warn > 0, ]
cat("files with fail, error or warn:", nrow(bad), "\n")
if (nrow(bad)) print(bad, row.names = FALSE)
nr <- grep("NO RESULT", x, value = TRUE)
cat("files without a RESULT line:", length(nr), "\n")
if (length(nr)) cat(nr, sep = "\n")
