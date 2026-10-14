# Lane illse: totals of a suite scan's .sum file (dev/illse-scan.sh),
# per package: files, pass, fail, error, skip, warn, and every file with
# a fail, error, warning or no RESULT line, with its DETAIL lines; and
# the lib: line of every per-file log (which frmtmb each loaded).
#   Rscript dev/illse-scansum.R <run>
run <- commandArgs(TRUE)[1]
d <- "C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-log/scan"
s <- readLines(file.path(d, paste0(run, ".sum")))
cat(tail(s, 1), "\n")
r <- s[grepl(" RESULT ", s)]
pkg <- sub(" .*", "", r)
num <- function(k) as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", r))
tb <- data.frame(pkg, pass = num("pass"), fail = num("fail"),
                 err = num("err"), skip = num("skip"), warn = num("warn"))
agg <- aggregate(cbind(files = 1, pass, fail, err, skip, warn) ~ pkg,
                 data = tb, FUN = sum)
print(agg, row.names = FALSE)
cat("total files", sum(agg$files), "pass", sum(agg$pass), "fail",
    sum(agg$fail), "err", sum(agg$err), "skip", sum(agg$skip), "warn",
    sum(agg$warn), "\n")
bad <- r[tb$fail > 0 | tb$err > 0 | tb$warn > 0]
nores <- s[grepl("NO RESULT", s)]
cat("\nfiles with a fail, error or warning:\n"); cat(bad, sep = "\n")
cat("\nno result:\n"); cat(nores, sep = "\n")
logs <- list.files(file.path(d, paste0(run, "-suite")), "[.]txt$",
                   full.names = TRUE)
libs <- vapply(logs, function(f) {
  l <- grep("^lib: ", readLines(f, warn = FALSE), value = TRUE)
  if (length(l)) sub(".*; frmtmb [^ ]+ ", "", l[1]) else "NONE"
}, "")
cat("\nfrmtmb loaded from:\n"); print(table(libs))
for (f in logs) {
  det <- grep("^DETAIL\t(FAIL|ERROR|WARN)", readLines(f, warn = FALSE),
              value = TRUE)
  if (length(det)) {
    cat("\n", basename(f), "\n", sep = "")
    cat(substr(det, 1, 300), sep = "\n")
  }
}
