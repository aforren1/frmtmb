# Per-package totals of a lane suite run, from its tier log
# (dev/ordmix-suite-<tier>/<tier>.log), with any rerun tier's files
# taking the place of their first result, and the lib: line of every
# per-file log checked. Usage: Rscript dev/ordmix-suite-sum.R <tier>
# [<rerun tier>]; prints the block the findings paste.
args <- commandArgs(TRUE)
root <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev"
rd <- function(tier) {
  x <- readLines(file.path(root, paste0("ordmix-suite-", tier),
                           paste0(tier, ".log")))
  ran <- grep("^RAN ", x, value = TRUE)
  miss <- grep("NO RESULT LINE", x, value = TRUE)
  x <- grep(" RESULT ", x, value = TRUE)
  m <- regmatches(x, regexec(paste0("^(\\S+) RESULT (\\S+) (\\S+) (\\S+) ",
                                    "pass=(\\d+) fail=(\\d+) err=(\\d+) ",
                                    "skip=(\\d+) warn=(\\d+)"), x))
  d <- do.call(rbind, lapply(m, function(v) {
    data.frame(pkg = v[2], file = v[5], pass = as.integer(v[6]),
               fail = as.integer(v[7]), err = as.integer(v[8]),
               skip = as.integer(v[9]), warn = as.integer(v[10]))
  }))
  list(d = d, ran = ran, miss = miss)
}
a <- rd(args[1])
cat(a$ran, "\n")
if (length(a$miss)) cat("missing:", a$miss, sep = "\n")
f <- a$d
if (length(args) >= 2) {
  r <- rd(args[2])
  cat("rerun", r$ran, "\n")
  key <- function(d) paste(d$pkg, d$file)
  f <- rbind(f[!key(f) %in% key(r$d), ], r$d)
}
tot <- aggregate(cbind(files, pass, fail, err, skip, warn) ~ pkg,
                 data = transform(f, files = 1L), FUN = sum)
print(tot, row.names = FALSE)
cat(sprintf("all: %d files, pass %d fail %d err %d skip %d warn %d\n",
            nrow(f), sum(f$pass), sum(f$fail), sum(f$err), sum(f$skip),
            sum(f$warn)))
cat("\nfiles with a failure, an error, a skip or a warning:\n")
print(f[f$fail + f$err + f$skip + f$warn > 0, ], row.names = FALSE)
# every per-file log must load frmtmb from the lane library
logs <- list.files(file.path(root, paste0("ordmix-suite-", args[1])),
                   pattern = "--.*[.]txt$", full.names = TRUE)
libs <- vapply(logs, function(l) {
  x <- grep("^lib:", readLines(l, n = 5, warn = FALSE), value = TRUE)
  if (length(x)) x[1] else NA_character_
}, "")
cat("\nlogs:", length(logs), " loading frmtmb from wt-ordmix-lib:",
    sum(grepl("frmtmb: C:/Users/adf44/source/r/wt-ordmix-lib/frmtmb",
              libs, fixed = TRUE)), "\n")
