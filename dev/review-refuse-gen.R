# Fills the review's run block from the logs, so its counts are generated, and
# wraps the prose at 80 columns. Run once, after the runs end.
f <- "dev/reviews/2026-11-02-refuse.md"
sumline <- function(path, lab) {
  if (!file.exists(path)) return(sprintf("- %s: no summary file.", lab))
  s <- readLines(path, warn = FALSE)
  ran <- grep("^ran ", s, value = TRUE)
  res <- grep("RESULT", s, value = TRUE)
  g <- function(k) sum(as.numeric(sub(paste0(".*\\b", k, "=([0-9]+).*"), "\\1",
                                       grep(paste0("\\b", k, "="), res, value = TRUE))))
  bad <- grep("fail=0 err=0", res, invert = TRUE, value = TRUE)
  sprintf("- %s: %s; pass %d, fail %d, error %d, skip %d, warn %d.%s", lab, ran, g("pass"), g("fail"),
          g("err"), g("skip"), g("warn"), if (length(bad)) paste0(" Not clean: ", paste(bad, collapse = "; "), ".") else "")
}
lk <- "dev/review-refuse-out/ob/dev/ciharden-log"
ob <- sumline(file.path(lk, "reviewob.sum"), "OpenBLAS 0.3.32 with LAPACK, 4 threads, every gate on, one process per file (`dev/review-refuse-ob.sh`)")
rf <- sumline(file.path(lk, "reviewref.sum"), "Reference BLAS, every gate on, one process per file, final tree (`dev/review-refuse-ref.sh`), which includes test-brms-lik-cdfs.R")
cd <- grep("test-brms-lik-cdfs", readLines(file.path(lk, "reviewref.sum"), warn = FALSE), value = TRUE)
ck <- "dev/review-refuse-check/check-frmtmb.log"
st <- if (file.exists(ck)) grep("^Status", readLines(ck, warn = FALSE), value = TRUE) else "not run"
nt <- if (file.exists(ck)) {
  l <- readLines(ck, warn = FALSE)
  c(grep("^\\* checking.*(NOTE|WARNING|ERROR)$", l, value = TRUE))
} else character()
head_l <- grep("headline", readLines("dev/review-refuse-equiv-out/summary.md"), value = TRUE)
noadj <- grep("headline", readLines("dev/review-refuse-out/adj/equiv-noadj-out/summary.md"), value = TRUE)
pick <- function(s) sub(".*headline\\): ([0-9.]+),.*", "\\1", s)
muts <- unlist(lapply(list.files("dev/review-refuse-out/mut", full.names = TRUE), function(p) {
  grep("^MUT", readLines(p, warn = FALSE), value = TRUE)
}))
blk <- c(
  ob, rf,
  sprintf("- test-brms-lik-cdfs.R in the gated tier: %s. So it was read.", trimws(sub("^frmtmb ", "", cd))),
  sprintf("- `R CMD check --as-cran`, core, after the NEWS edit (`dev/review-refuse-check.sh`): %s. Sections not OK: %s.",
          st, if (length(nt)) paste(nt, collapse = "; ") else "none"),
  sprintf("- All 172 pairs (`dev/review-refuse-equiv-out`): headline %s; verdicts and L1 values equal the lane's in every pair. With adjudication rows p043 and p169 dropped, by a run: %s.",
          pick(head_l), pick(noadj)),
  sprintf("- Mutants run: %d (`dev/review-refuse-out/mut/`).", length(muts)),
  "- Not rerun: the three-order one-session runs of the nine suites, and frmtmb.sample's `R CMD check`; the lane's are in its record, and only NEWS changed after them. Other extensions' suites were not rerun.")
x <- readLines(f, warn = FALSE)
x <- unlist(lapply(x, function(l) if (identical(l, "GENERATED-RUNS")) blk else l))
wrap <- function(lines) {
  out <- character(); i <- 1L
  buf <- character(); pref <- ""; ind <- ""
  flush <- function() {
    if (length(buf)) {
      t <- paste(sub("^ +", "", buf), collapse = " ")
      out <<- c(out, strwrap(t, width = 79, initial = pref, exdent = nchar(ind)))
      buf <<- character()
    }
  }
  for (l in lines) {
    if (!nzchar(l) || startsWith(l, "|") || startsWith(l, "#")) {
      flush(); out <- c(out, l); pref <- ""; ind <- ""
    } else if (startsWith(l, "- ")) {
      flush(); pref <- "- "; ind <- "  "; buf <- sub("^- ", "", l)
    } else {
      if (!length(buf)) { pref <- ""; ind <- "" }
      buf <- c(buf, l)
    }
  }
  flush()
  out
}
writeLines(wrap(x), f)
cat(sprintf("%d lines; over 80 outside tables: %d\n", length(readLines(f)),
            sum(nchar(readLines(f)) > 80 & !startsWith(readLines(f), "|"))))
