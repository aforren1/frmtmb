# Lane wt-conditions (item 2.6e): convert every stop(), warning() and
# message() that frmtmb and its extensions raise to frm_stop(),
# frm_warning() and frm_message(), driven by the parse tree.
#
#   Rscript dev/conditions-rewrite.R rewrite <snapshot-dir>
#       copies each R/ tree to <snapshot-dir> first, then rewrites in
#       place and prints the counts
#   Rscript dev/conditions-rewrite.R verify <snapshot-dir>
#       proves that every file differs from its snapshot ONLY by the
#       transformation below: the expressions of the snapshot, put
#       through to_new(), deparse identically to the expressions of the
#       current file. Line breaks and indentation are free to differ,
#       which is what lets a line pushed past 80 columns be reflowed by
#       hand without the proof losing its grip.
#
# Run from the worktree root. The transformation, decided site by site
# in dev/conditions-findings.md:
#   * stop / warning / message  ->  frm_stop / frm_warning / frm_message,
#     head renamed, arguments untouched;
#   * warning(warningCondition(<args>))  ->
#     frm_warning(<args>, call. = FALSE), because warningCondition()
#     records no call and frm_warning() takes `class` itself;
#   * stop(structure(class = c("frmtmb_fit_error", "simpleError",
#     "error", "condition"), list(message = msg, call = NULL)))  ->
#     frm_stop(msg, call. = FALSE, class = "frmtmb_fit_error");
#   * stop(e), a rethrow of a caught condition, is left alone: it keeps
#     whatever class the condition already has.
# R/conditions.R, which defines the helpers, is not rewritten.

av <- commandArgs(trailingOnly = TRUE)
mode <- av[1L]
snap <- av[2L]
stopifnot(mode %in% c("rewrite", "verify", "restore"), !is.na(snap))

roots <- c("R", file.path(sort(list.dirs("extensions", recursive = FALSE)),
                          "R"))
skip_files <- "R/conditions.R"
heads <- c(stop = "frm_stop", warning = "frm_warning",
           message = "frm_message")

fit_error_old <- quote(stop(structure(
  class = c("frmtmb_fit_error", "simpleError", "error", "condition"),
  list(message = msg, call = NULL))))
fit_error_new <- quote(frm_stop(msg, call. = FALSE,
                                class = "frmtmb_fit_error"))

is_rethrow <- function(e) {
  is.call(e) && identical(e[[1L]], as.name("stop")) && length(e) == 2L &&
    is.name(e[[2L]])
}
is_wc <- function(e) {
  is.call(e) && identical(e[[1L]], as.name("warning")) && length(e) == 2L &&
    is.call(e[[2L]]) && identical(e[[2L]][[1L]], as.name("warningCondition"))
}

# the AST transformation, one place, used by the proof
to_new <- function(e) {
  if (!is.call(e)) return(e)
  if (identical(e, fit_error_old)) return(fit_error_new)
  if (is_rethrow(e)) return(e)
  if (is_wc(e)) {
    inner <- as.list(e[[2L]])[-1L]
    e <- as.call(c(list(as.name("frm_warning")), inner,
                   list(call. = FALSE)))
  } else if (is.name(e[[1L]]) && as.character(e[[1L]]) %in% names(heads)) {
    e[[1L]] <- as.name(heads[[as.character(e[[1L]])]])
  }
  for (i in seq_along(e)) {
    # an empty argument (`x[, 1]`) errors when touched
    a <- tryCatch(if (is.call(e[[i]])) e[[i]], error = function(err) NULL)
    if (!is.null(a)) e[[i]] <- to_new(a)
  }
  e
}

files <- unlist(lapply(roots, function(r) {
  list.files(r, pattern = "[.][Rr]$", full.names = TRUE)
}))
# the worked-example families shipped in inst/ (sourced by the vignettes
# and by the tests) raise refusals the suite asserts, so they are
# frmtmb's conditions too; the sweep found them, in a second pass
inst_files <- unlist(lapply(c("inst", file.path(sort(list.dirs(
  "extensions", recursive = FALSE)), "inst")), function(r) {
    if (dir.exists(r)) list.files(r, pattern = "[.][Rr]$", full.names = TRUE,
                                  recursive = TRUE)
  }))
files <- setdiff(c(files, inst_files), skip_files)

if (mode == "restore") {
  # put the snapshot back and remove it, so a rewrite can run again
  for (f in files) stopifnot(file.copy(file.path(snap, f), f,
                                       overwrite = TRUE))
  unlink(snap, recursive = TRUE)
  cat("restored", length(files), "files\n")
  quit(status = 0L)
}

if (mode == "rewrite") {
  # a file with a snapshot was rewritten by an earlier pass and is left
  # alone, so a pass that adds files does not rewrite a file twice
  files_all <- files
  files <- files[!file.exists(file.path(snap, files))]
  cat("rewriting", length(files), "files not yet snapshotted\n")
  for (f in files) {
    dir.create(dirname(file.path(snap, f)), recursive = TRUE,
               showWarnings = FALSE)
    stopifnot(file.copy(f, file.path(snap, f), overwrite = FALSE))
  }
  counts <- list()
  bump <- function(pkg, what) {
    key <- paste(pkg, what)
    counts[[key]] <<- (if (is.null(counts[[key]])) 0L else
      counts[[key]]) + 1L
  }
  for (f in files) {
    pkg <- if (startsWith(f, "R/")) "frmtmb" else
      if (startsWith(f, "inst/")) "frmtmb-inst" else
        sub("^extensions/([^/]+)/(R|inst)/.*$", "\\1", f)
    src <- readLines(f, warn = FALSE, encoding = "UTF-8")
    pd <- getParseData(parse(f, keep.source = TRUE, encoding = "UTF-8"),
                       includeText = TRUE)
    edits <- list()   # line, col1, line2, col2, text
    fc <- pd[pd$token == "SYMBOL_FUNCTION_CALL" &
               pd$text %in% names(heads), ]
    for (k in seq_len(nrow(fc))) {
      tok <- fc[k, ]
      head_expr <- pd[pd$id == tok$parent, ]
      call_expr <- pd[pd$id == head_expr$parent, ]
      e <- str2lang(call_expr$text)
      if (is_rethrow(e)) {
        bump(pkg, "exempt stop(<condition>) rethrow")
        cat("EXEMPT rethrow ", f, ":", tok$line1, " ", call_expr$text,
            "\n", sep = "")
        next
      }
      stopifnot(identical(substr(src[tok$line1], tok$col1, tok$col2),
                          tok$text))
      if (identical(e, fit_error_old)) {
        edits[[length(edits) + 1L]] <- list(special = TRUE,
          l1 = call_expr$line1, c1 = call_expr$col1,
          l2 = call_expr$line2, c2 = call_expr$col2,
          text = "frm_stop(msg, call. = FALSE, class = \"frmtmb_fit_error\")")
        bump(pkg, "stop(structure(frmtmb_fit_error)) special")
        next
      }
      if (is_wc(e)) {
        # the text between warningCondition( and its closing paren is
        # kept byte for byte, so the message arguments are not re-typed
        kids <- pd[pd$parent == call_expr$id, ]
        wc <- kids[kids$token == "expr" & kids$line1 * 1e4 + kids$col1 >
                     tok$line1 * 1e4 + tok$col2, ][1L, ]
        wk <- pd[pd$parent == wc$id, ]
        op <- wk[wk$token == "'('", ][1L, ]
        cp <- wk[wk$token == "')'", ]
        cp <- cp[nrow(cp), ]
        body <- if (op$line1 == cp$line1) {
          substr(src[op$line1], op$col1 + 1L, cp$col1 - 1L)
        } else {
          paste(c(substring(src[op$line1], op$col1 + 1L),
                  if (cp$line1 - op$line1 > 1L)
                    src[(op$line1 + 1L):(cp$line1 - 1L)],
                  substr(src[cp$line1], 1L, cp$col1 - 1L)),
                collapse = "\n")
        }
        edits[[length(edits) + 1L]] <- list(special = TRUE,
          l1 = call_expr$line1, c1 = call_expr$col1,
          l2 = call_expr$line2, c2 = call_expr$col2,
          text = paste0("frm_warning(", body, ", call. = FALSE)"))
        bump(pkg, "warning(warningCondition()) special")
        next
      }
      # a plain site: rename the head token only
      edits[[length(edits) + 1L]] <- list(special = FALSE,
        l1 = tok$line1, c1 = tok$col1, l2 = tok$line1, c2 = tok$col2,
        text = heads[[tok$text]])
      # continuation lines aligned under the arguments move right with
      # the longer name, so the call keeps its alignment; a string
      # literal spanning lines would change here, and the proof below
      # compares string contents, so it would catch that
      grow <- nchar(heads[[tok$text]]) - nchar(tok$text)
      if (call_expr$line2 > tok$line1) {
        for (ln in (tok$line1 + 1L):call_expr$line2) {
          ind <- nchar(sub("^( *).*$", "\\1", src[ln]))
          if (nzchar(trimws(src[ln])) && ind >= tok$col2 + 1L) {
            edits[[length(edits) + 1L]] <- list(special = FALSE,
              l1 = ln, c1 = 1L, l2 = ln, c2 = 0L,
              text = strrep(" ", grow))
          }
        }
      }
      bump(pkg, tok$text)
    }
    if (!length(edits)) next
    # a special edit spans a whole call; drop the head renames nested in
    # it, which its own text already carries
    spans <- Filter(function(x) x$special, edits)
    inside <- function(x) {
      any(vapply(spans, function(s) {
        !x$special &&
          (x$l1 > s$l1 || (x$l1 == s$l1 && x$c1 >= s$c1)) &&
          (x$l2 < s$l2 || (x$l2 == s$l2 && x$c2 <= s$c2))
      }, TRUE))
    }
    edits <- Filter(Negate(inside), edits)
    ord <- order(-vapply(edits, `[[`, 1, "l1"),
                 -vapply(edits, `[[`, 1, "c1"))
    for (x in edits[ord]) {
      pre <- substr(src[x$l1], 1L, x$c1 - 1L)
      post <- substring(src[x$l2], x$c2 + 1L)
      new <- strsplit(paste0(pre, x$text, post), "\n", fixed = TRUE)[[1L]]
      src <- c(src[seq_len(x$l1 - 1L)], new,
               src[seq_along(src) > x$l2])
    }
    con <- file(f, open = "wb")
    writeLines(enc2utf8(src), con, sep = "\n", useBytes = TRUE)
    close(con)
  }
  cat("\n== counts by package ==\n")
  for (k in sort(names(counts))) cat(sprintf("%-60s %5d\n", k, counts[[k]]))
}

if (mode == "rewrite") files <- files_all

# the proof, run after a rewrite and again after any hand reflow
bad <- 0L
long <- 0L
for (f in files) {
  old <- parse(file.path(snap, f), keep.source = FALSE, encoding = "UTF-8")
  new <- parse(f, keep.source = FALSE, encoding = "UTF-8")
  same <- length(old) == length(new) &&
    all(vapply(seq_along(old), function(i) {
      identical(deparse(to_new(old[[i]])), deparse(new[[i]]))
    }, TRUE))
  if (!same) {
    bad <- bad + 1L
    cat("MISMATCH ", f, "\n", sep = "")
  }
  # only lines the rename pushed past 80 columns count: some rule tables
  # were over 80 before this lane and are not its business
  wid <- function(p) {
    x <- readLines(p, warn = FALSE, encoding = "UTF-8")
    sort(x[nchar(x, type = "width") > 80L])
  }
  was <- trimws(wid(file.path(snap, f)))
  now <- wid(f)
  now_renamed <- gsub("frm_(stop|warning|message)[(]", "\\1(", trimws(now))
  for (i in which(!now_renamed %in% was)) {
    long <- long + 1L
    cat("OVER80 ", f, ": ", substr(trimws(now[i]), 1, 60), "\n", sep = "")
  }
}
cat("files compared ", length(files), ", mismatched ", bad,
    ", lines over 80 columns ", long, "\n", sep = "")
if (bad) quit(status = 1L)
