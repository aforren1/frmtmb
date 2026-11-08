# Summarize dev/nlstart-scan.R output. The reference of a row is the best
# log-likelihood any arm reached on it, by default or by random start.
# A row is short when its default fit is below the reference by more
# than 1e-4; silent when the default fit raised no warning of any class.
# Usage: Rscript dev/nlstart-scansum.R <dir> [<dir> ...]
dirs <- commandArgs(trailingOnly = TRUE)
read_arm <- function(dir) {
  fs <- list.files(dir, "[.]txt$", full.names = TRUE)
  l <- unlist(lapply(fs, function(f) grep("^ROW ", readLines(f), value = TRUE)))
  x <- read.table(text = sub("^ROW ", "", l), stringsAsFactors = FALSE,
                  col.names = c("family", "seed", "slope", "body", "default",
                                "ref", "conv", "t_def", "nwarn", "wclass",
                                "nmsg"))
  x$key <- paste(x$family, x$seed, x$slope, x$body)
  x
}
arms <- lapply(dirs, read_arm)
names(arms) <- basename(dirs)
keys <- Reduce(intersect, lapply(arms, `[[`, "key"))
best <- sapply(keys, function(k) {
  max(unlist(lapply(arms, function(a) {
    r <- a[a$key == k, ]
    c(r$default, r$ref)
  })), na.rm = TRUE)
})
cat("rows in every arm:", length(keys), "\n")
for (nm in names(arms)) {
  a <- arms[[nm]]
  a <- a[match(keys, a$key), ]
  a$short <- best[a$key] - a$default
  short <- !is.na(a$short) & a$short > 1e-4
  na <- is.na(a$default)
  silent <- a$nwarn == 0
  cat(sprintf(paste0("%s: rows %d, default NA %d, short > 1e-4: %d ",
                     "(silent %d, warned %d), short > 1: %d (silent %d), ",
                     "wall time of default fits %.1f s\n"),
              nm, nrow(a), sum(na), sum(short), sum(short & silent),
              sum(short & !silent), sum(short & a$short > 1),
              sum(short & a$short > 1 & silent), sum(a$t_def)))
  o <- a[short, c("family", "seed", "slope", "body", "default", "conv",
                  "nwarn", "wclass", "nmsg")]
  o$best <- best[a$key[short]]
  o$short <- round(a$short[short], 4)
  print(o[order(-o$short), ], row.names = FALSE)
}
if (length(arms) == 2L) {
  a <- arms[[1]][match(keys, arms[[1]]$key), ]
  b <- arms[[2]][match(keys, arms[[2]]$key), ]
  d <- b$default - a$default
  cat("\n", names(arms)[2], "minus", names(arms)[1], "default logLik:",
      "higher by > 1e-4:", sum(d > 1e-4, na.rm = TRUE),
      " lower by > 1e-4:", sum(d < -1e-4, na.rm = TRUE),
      " within 1e-4:", sum(abs(d) <= 1e-4, na.rm = TRUE), "\n")
  lo <- which(d < -1e-4)
  if (length(lo)) {
    print(data.frame(key = keys[lo], first = a$default[lo],
                     second = b$default[lo], best = best[keys[lo]]),
          row.names = FALSE)
  }
  cat("time ratio", names(arms)[2], "/", names(arms)[1], ":",
      round(sum(b$t_def) / sum(a$t_def), 3), "\n")
  fam_t <- tapply(b$t_def, b$family, sum) / tapply(a$t_def, a$family, sum)
  print(round(fam_t, 3))
}
