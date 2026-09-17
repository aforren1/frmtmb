# Lane wt-conditions: the generated blocks of dev/conditions-findings.md,
# read from the run directories. Nothing here is typed from memory.
#   Rscript dev/conditions-sweep-summary.R <rundir>... > <out>
dirs <- commandArgs(trailingOnly = TRUE)

num_after <- function(x, k) {
  as.integer(vapply(regmatches(x, regexec(paste0(k, " +([0-9]+)"), x)),
                    function(m) if (length(m) > 1L) m[2L] else NA_character_,
                    ""))
}
cat("== suites, one file per R process ==\n")
cat("runner dev/conditions-sweep-onefile.R, sums failures AND errors\n")
for (d in dirs) {
  s <- readLines(file.path(d, "summary.txt"), warn = FALSE)
  pkg <- ifelse(grepl("extensions/", s),
                sub("^.*extensions/([^/]+)/.*$", "\\1", s), "frmtmb")
  for (p in sort(unique(pkg))) {
    x <- s[pkg == p]
    ok <- grepl("^BLOCKS", x)
    cat(sprintf(paste0("%-16s files %3d (no result line %d) blocks %4d ",
                       "PASS %5d FAIL %d ERROR %d SKIP %d\n"),
                p, length(x), sum(!ok), sum(num_after(x[ok], "BLOCKS")),
                sum(num_after(x[ok], "PASS")), sum(num_after(x[ok], "FAIL")),
                sum(num_after(x[ok], "ERROR")),
                sum(num_after(x[ok], "SKIP"))))
  }
  bad <- s[!grepl("^BLOCKS", s) | num_after(s, "FAIL") > 0 |
             num_after(s, "ERROR") > 0]
  if (length(bad)) cat("  not clean:\n", paste0("    ", bad, "\n"), sep = "")
  for (st in c(".start", ".end")) {
    f <- file.path(d, st)
    if (file.exists(f)) cat("  ", st, readLines(f), "\n")
  }
}

# The literal templates of every frm_stop() call, to tell an unclassed
# error frmtmb raised (a hole in the conversion) from one it did not
roots <- c("R", file.path(list.dirs("extensions", recursive = FALSE), "R"))
lits <- character(0)
walk <- function(e) {
  if (!is.call(e)) return(invisible())
  if (is.name(e[[1L]]) &&
      as.character(e[[1L]]) %in% c("frm_stop", "stop", "stopifnot")) {
    for (i in seq_along(e)[-1L]) {
      v <- tryCatch(if (is.character(e[[i]])) e[[i]], error = function(x) NULL)
      if (length(v) == 1L && nchar(v) >= 20L) lits <<- c(lits, v)
    }
  }
  for (i in seq_along(e)) {
    a <- tryCatch(if (is.call(e[[i]])) e[[i]], error = function(x) NULL)
    if (!is.null(a)) walk(a)
  }
}
for (r in roots) for (f in list.files(r, "[.][Rr]$", full.names = TRUE)) {
  for (e in parse(f, keep.source = FALSE)) walk(e)
}
lits <- unique(lits)

cat("\n== sweep: conditions the suites' own expect_error() caught ==\n")
rows <- do.call(rbind, lapply(dirs, function(d) {
  fs <- list.files(d, pattern = "[.]tsv$", full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) {
    utils::read.delim(f, stringsAsFactors = FALSE, quote = "\"")
  }))
}))
rows$pkg <- ifelse(grepl("extensions/", rows$file),
                   sub("^.*extensions/([^/]+)/.*$", "\\1", rows$file),
                   "frmtmb")
rows$frmtmb <- grepl("(^|/)frmtmb_error(/|$)", rows$class)
cat(sprintf("%-16s %6s %13s %6s\n", "package", "caught", "frmtmb_error",
            "other"))
for (p in c(sort(unique(rows$pkg)), "ALL")) {
  x <- if (p == "ALL") rows else rows[rows$pkg == p, ]
  cat(sprintf("%-16s %6d %13d %6d\n", p, nrow(x), sum(x$frmtmb),
              sum(!x$frmtmb)))
}
own <- rows[rows$frmtmb, ]
own_sub <- sub("/frmtmb_error.*$", "", own$class)
cat("\nfrmtmb_error class vectors seen:\n")
print(table(own$class))

other <- rows[!rows$frmtmb, ]
hit <- vapply(other$message, function(m) {
  any(vapply(lits, function(l) grepl(l, m, fixed = TRUE), TRUE))
}, TRUE)
other$template <- hit
cat("\nnot frmtmb_error, by class and whether the message carries a",
    "literal\nof a stop()/frm_stop()/stopifnot() call in frmtmb's",
    "sources:\n")
print(table(class = other$class, frmtmb_literal = other$template))
cat("\nevery distinct non-frmtmb_error message (class | package | text):\n")
key <- unique(other[, c("class", "pkg", "message", "template")])
key <- key[order(key$template, key$class, key$pkg), ]
for (i in seq_len(nrow(key))) {
  cat(sprintf("[%s] %s | %s | %s\n",
              if (key$template[i]) "LITERAL" else "foreign", key$class[i],
              key$pkg[i], substr(key$message[i], 1, 150)))
}

cat("\nnot frmtmb_error, by origin (patterns below, first match wins):\n")
origin <- function(cls, msg) {
  if (grepl("brms_error", cls)) return("brms itself (the agreement tests)")
  pats <- c(
    "should be one of" = "base match.arg() inside a frmtmb function",
    "error in evaluating the argument" =
      "R's S4 dispatch rewraps a frmtmb_error as simpleError",
    "unused argument|not found$|object '.*' not found" =
      "R evaluation: an argument or object that does not exist",
    "invalid type|variable lengths differ|has new levels" =
      "stats::model.frame() and its relatives",
    "link not recognised" = "stats::make.link()",
    "missing in the draws object" = "the posterior package",
    "^bad x: 3$|^brms parameterizes a|translation rule" =
      "test code: a test helper or a deliberate absent case")
  for (p in names(pats)) if (grepl(p, msg)) return(pats[[p]])
  "UNCLASSIFIED"
}
other$origin <- mapply(origin, other$class, other$message)
tab <- sort(table(other$origin), decreasing = TRUE)
for (k in names(tab)) cat(sprintf("  %3d  %s\n", tab[[k]], k))
