# Summary of the warning-precedence log of dev/rel075-warnscan.sh
# (dev/rel070-profile.R): per fit, the parameters each report names,
# and every parameter two reports name, or one report names twice.
#
#   Rscript dev/rel075-prof-sum.R > dev/rel075-log/prof-sum.txt
fs <- list.files(c("dev/rel075-log/prof-ext", "dev/rel075-log/prof-core"),
                 pattern = "^[0-9]+[.]txt$",
                 full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) {
  l <- readLines(f, warn = FALSE)
  l <- l[startsWith(l, "REP\t")]
  if (!length(l)) return(NULL)
  s <- strsplit(l, "\t", fixed = TRUE)
  data.frame(file = vapply(s, `[`, "", 2), pid = vapply(s, `[`, "", 3),
             id = vapply(s, `[`, "", 4), kind = vapply(s, `[`, "", 5),
             pars = vapply(s, function(v) if (length(v) >= 6) v[6] else "",
                           ""))
}))
# keyed by file too: Windows reuses a PID, and two test files run by
# processes with one PID append to one log (0.70.0 run: pid 11132)
x$fit <- paste(x$file, x$pid, x$id)
cat("processes with a log:", length(fs), "\n")
cat("files with a log line:", length(unique(x$file)), "\n")
cat("fits (fit_assembled() calls):", sum(x$kind == "FIT"), "\n")
cat("lines with no fit id (outside any fit):", sum(x$id == "0"), "\n")
tab <- table(x$kind[x$kind != "FIT"])
for (k in names(tab)) cat(sprintf("  %-15s %d\n", k, tab[[k]]))
rep_kinds <- c("SE", "BOUNDARY", "SEPARATION", "NONFINITE", "FLAT")
r <- x[x$kind %in% rep_kinds, ]
cat("fits with at least one report:", length(unique(r$fit)), "\n")
nk <- tapply(r$kind, r$fit, function(k) length(unique(k)))
cat("fits with two or more kinds of report:", sum(nk >= 2), "\n")
bad <- 0L
for (fi in unique(r$fit)) {
  y <- r[r$fit == fi, ]
  # one row per (kind, parameter); FLAT names no parameter, and stands
  # for every parameter but a bound-held one (se_report()'s rule)
  pp <- do.call(rbind, lapply(seq_len(nrow(y)), function(i) {
    p <- strsplit(y$pars[i], ",", fixed = TRUE)[[1]]
    if (!length(p)) p <- "<all>"
    data.frame(kind = y$kind[i], par = sub(":.*", "", p),
               why = ifelse(grepl(":", p), sub(".*:", "", p), ""))
  }))
  dup <- pp$par[duplicated(pp$par)]
  flat_se <- any(pp$kind == "FLAT") &&
    any(pp$kind == "SE" & pp$why != "bound")
  if (length(dup) || flat_se) {
    bad <- bad + 1L
    cat("TWICE", y$file[1], "fit", fi, ":",
        paste(paste0(pp$kind, "=", pp$par, ifelse(nzchar(pp$why),
                                                  paste0("(", pp$why, ")"),
                                                  "")), collapse = " "),
        "\n")
  }
  if (nk[[fi]] >= 2) {
    cat("MULTI", y$file[1], "fit", fi, ":",
        paste(paste0(pp$kind, "=", pp$par), collapse = " "), "\n")
  }
}
cat("fits where a parameter got two reports:", bad, "\n")
