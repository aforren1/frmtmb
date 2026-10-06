# Punch round 1, B1: per model, fits, errors and how many fits raised
# each kind of warning, for a dev/ordmix-*-fa-log directory.
# Usage: Rscript dev/ordmix-p1-fa-sum.R <log dir>
dir <- commandArgs(TRUE)[1]
kinds <- c(flat = paste0("flat|not identified|no usable standard error|",
                         "one distribution|not in the likelihood"),
           # the log cuts warnings at 200 characters, so the
           # degenerate warning is known by its opening
           degenerate = "[)]: component [0-9]+ [(]",
           gradient = "gradient at the reported optimum is not finite",
           other = "")
rows <- NULL
for (f in list.files(dir, "[.]txt$", full.names = TRUE)) {
  l <- readLines(f, warn = FALSE)
  rep <- grep("^REP ", l)
  fits <- 0L; err <- 0L; hit <- setNames(integer(4), names(kinds))
  for (i in seq_along(rep)) {
    if (grepl(" ERROR ", l[rep[i]])) { err <- err + 1L; next }
    fits <- fits + 1L
    end <- if (i < length(rep)) rep[i + 1L] - 1L else length(l)
    w <- grep("^   WARN", l[seq_len(end - rep[i]) + rep[i]], value = TRUE)
    w <- w[!grepl("[se]", w, fixed = TRUE)]
    seen <- rep(FALSE, length(w))
    for (k in names(kinds)[1:3]) {
      m <- grepl(kinds[[k]], w)
      hit[[k]] <- hit[[k]] + any(m)
      seen <- seen | m
    }
    hit[["other"]] <- hit[["other"]] + any(!seen)
  }
  rows <- rbind(rows, data.frame(model = sub("[.]txt$", "", basename(f)),
                                 fits = fits, errors = err,
                                 as.list(hit)))
}
print(rows, row.names = FALSE)
