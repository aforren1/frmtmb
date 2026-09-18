# False-alarm rate of the three refusals this lane adds, with brms
# 2.23.0 as the judge.
#
#   Rscript dev/adefects-falsealarm.R > dev/adefects-log/falsealarm.txt 2>&1
#
# WHAT IS HARVESTED. Every `ar()`, `ma()`, `arma()`, `cosy()` and
# `unstr()` call and every `hypothesis()` string that appears in:
#
#   * brms 2.23.0's own tests, vignettes and Rd examples
#     (dev/brms-suite/brms, re-fetched per dev/brms-suite-audit.md
#     section 1),
#   * this repository's test suites, vignettes and Rd examples.
#
# The harvest PARSES each file and walks the syntax tree, so a term
# inside a comment or a string is not counted and a call split over
# lines is counted once. A grep for `ar(` scored 49 copies of
# `get_dpar(prep, "mu", i = i)` on its first try, which is why this
# parses.
#
# WHAT IS ASKED. Each harvested call goes to brms's own constructor
# (`brms::ar()` and friends, which validate and do not evaluate their
# arguments) and to frmtmb's parser
# (`frmtmb:::parse_autocor_call()`), and each harvested hypothesis
# string goes to `brms::hypothesis()` on brms's own stored example fit
# and to `frmtmb:::hyp_parse()`. Outcomes:
#
#   both accept        - agreement
#   both refuse        - agreement
#   FALSE ALARM        - frmtmb refuses what brms accepts
#   MISSED             - frmtmb accepts what brms refuses
#
# The brms grammar refusals are told apart from brms's other errors by
# their text, so a hypothesis brms rejects because the PARAMETER is
# unknown counts as grammar-accepted, which is what is being measured.
#
# ARM. The same harvest runs twice: once with the BASE build and the
# base tree (ADEFECTS_LIB = rellib-r3, ADEFECTS_ROOT = a `git archive`
# extract of f8b45ef) and once with this lane's. Diffing the two TSVs is
# what attributes a refusal to THIS LANE rather than to a divergence the
# package already had, such as the documented `cov = FALSE` refusal,
# which otherwise swamps the count.
lib <- Sys.getenv("ADEFECTS_LIB", "C:/Users/adf44/source/r/adefects-lib")
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
tag <- Sys.getenv("ADEFECTS_TAG", "lane")
root <- Sys.getenv("ADEFECTS_ROOT", ".")
cat("arm", tag, "lib", lib, "root", root, "\n")
cat("frmtmb", format(packageVersion("frmtmb")), "brms",
    format(packageVersion("brms")), "\n")

# brms's source is gitignored, so it is always read from this worktree
# whichever tree is being harvested; its sha256 is the recorded one.
brms_src <- "dev/brms-suite/brms"
stopifnot(dir.exists(brms_src))

# ------------------------------------------------------------- harvest

# R code of one file: an .R verbatim, an .Rmd's chunks, an .Rd's
# \examples{} section.
code_of <- function(f) {
  if (grepl("[.]R$", f)) return(parse(f, keep.source = FALSE))
  if (grepl("[.]Rmd$", f)) {
    ln <- readLines(f, warn = FALSE)
    open <- grep("^[`]{3}[{]r", ln)
    close <- grep("^[`]{3}\\s*$", ln)
    out <- character(0)
    for (o in open) {
      cl <- close[close > o]
      if (!length(cl)) next
      out <- c(out, ln[(o + 1L):(cl[1L] - 1L)])
    }
    return(parse(text = out, keep.source = FALSE))
  }
  rd <- tools::parse_Rd(f)
  tg <- vapply(rd, function(x) attr(x, "Rd_tag") %||% "", "")
  ex <- rd[tg == "\\examples"]
  if (!length(ex)) return(expression())
  txt <- paste(unlist(lapply(ex, function(e) as.character(unlist(e)))),
               collapse = "")
  parse(text = txt, keep.source = FALSE)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

files_of <- function(dirs, pats) {
  out <- character(0)
  for (d in dirs) {
    if (!dir.exists(d)) next
    for (p in pats) {
      out <- c(out, list.files(d, p, full.names = TRUE, recursive = TRUE))
    }
  }
  unique(out)
}

brms_files <- files_of(
  c(file.path(brms_src, "tests"), file.path(brms_src, "vignettes"),
    file.path(brms_src, "man")),
  c("[.]R$", "[.]Rmd$", "[.]Rd$"))
own_files <- files_of(
  c(file.path(root, "tests"), file.path(root, "vignettes"),
    file.path(root, "man"), file.path(root, "extensions")),
  c("[.]R$", "[.]Rmd$", "[.]Rd$"))
own_files <- own_files[!grepl("brms-suite", own_files, fixed = TRUE)]
# THIS LANE'S OWN TEST FILE IS NOT A DESIGN THE FIELD WRITES. It is
# built to make each new refusal fire, so counting it would report the
# lane's own constructions as false alarms and inflate the corpus. It is
# excluded from both arms; the base tree does not contain it anyway, and
# without this the two corpora differ by eight autocor calls.
own_files <- own_files[!grepl("test-adefects[.]R$", own_files)]

ac_fns <- c("ar", "ma", "arma", "cosy", "unstr")
ac_hits <- list()
hyp_hits <- list()

walk <- function(e, src) {
  if (!is.call(e)) return(invisible(NULL))
  fn <- e[[1L]]
  nm <- if (is.name(fn)) {
    as.character(fn)
  } else if (is.call(fn) && length(fn) == 3L &&
               as.character(fn[[1L]]) %in% c("::", ":::")) {
    as.character(fn[[3L]])
  } else ""
  if (nm %in% ac_fns) {
    ac_hits[[length(ac_hits) + 1L]] <<- list(call = e, src = src, fn = nm)
  }
  if (nm == "hypothesis" && length(e) >= 3L) {
    a <- e[[3L]]
    lit <- if (is.character(a)) a else if (
      is.call(a) && identical(a[[1L]], as.name("c")) &&
        all(vapply(as.list(a)[-1L], is.character, NA))) {
      unlist(as.list(a)[-1L])
    } else NULL
    if (!is.null(lit)) {
      hyp_hits[[length(hyp_hits) + 1L]] <<- list(h = lit, src = src)
    }
  }
  for (i in seq_along(e)) {
    p <- tryCatch({
      x <- e[[i]]
      if (is.call(x)) x else NULL
    }, error = function(err) NULL)
    if (!is.null(p)) walk(p, src)
  }
  invisible(NULL)
}

harvest <- function(fs, tag) {
  for (f in fs) {
    ex <- tryCatch(code_of(f), error = function(e) NULL)
    if (is.null(ex)) next
    for (i in seq_along(ex)) walk(ex[[i]], paste0(tag, ":", basename(f)))
  }
}
harvest(brms_files, "brms")
harvest(own_files, "own")

# An autocor call is a TERM only when it is a formula term: brms and
# frmtmb both read it from a formula, and `ar(x)` inside a plain
# expression is stats::ar(). Every harvested call is kept anyway, and
# the classification below does the work; a call whose function is one
# of the five names but whose arguments are not brms's grammar (say
# `ma(ns, 4)`, matrix arithmetic) shows up as an agreement, because
# brms refuses it too.
cat("harvested autocor-shaped calls:", length(ac_hits), "\n")
cat("harvested hypothesis strings:",
    sum(vapply(hyp_hits, function(x) length(x$h), 1L)), "\n")

# ------------------------------------------------------- the two judges

brms_ac <- function(cl) {
  cl2 <- cl
  cl2[[1L]] <- call("::", as.name("brms"), as.name(as.character(cl[[1L]])[
    length(as.character(cl[[1L]]))]))
  tryCatch({
    eval(cl2)
    TRUE
  }, error = function(e) structure(FALSE, msg = conditionMessage(e)))
}
frm_ac <- function(cl) {
  cl2 <- cl
  if (is.call(cl2[[1L]])) cl2[[1L]] <- as.name(as.character(cl[[1L]])[
    length(as.character(cl[[1L]]))])
  tryCatch({
    frmtmb:::parse_autocor_call(cl2, new.env())
    TRUE
  }, error = function(e) structure(FALSE, msg = conditionMessage(e)))
}

tab <- data.frame(src = character(0), call = character(0),
                  brms = logical(0), frm = logical(0))
for (h in ac_hits) {
  b <- brms_ac(h$call)
  f <- frm_ac(h$call)
  tab <- rbind(tab, data.frame(src = h$src, call = deparse1(h$call),
                               brms = isTRUE(b), frm = isTRUE(f)))
}
tab <- tab[!duplicated(tab$call), ]

cat("\n== D1/D2 autocor term grammar, distinct calls:", nrow(tab), "\n")
cat("  both accept :", sum(tab$brms & tab$frm), "\n")
cat("  both refuse :", sum(!tab$brms & !tab$frm), "\n")
cat("  FALSE ALARM (frmtmb refuses, brms accepts):",
    sum(tab$brms & !tab$frm), "\n")
for (i in which(tab$brms & !tab$frm)) {
  cat("    ", tab$src[i], " ", tab$call[i], "\n", sep = "")
}
cat("  MISSED (frmtmb accepts, brms refuses):",
    sum(!tab$brms & tab$frm), "\n")
for (i in which(!tab$brms & tab$frm)) {
  cat("    ", tab$src[i], " ", tab$call[i], "\n", sep = "")
}
cat("  accepted by both, listed:\n")
for (i in which(tab$brms & tab$frm)) {
  cat("    ", tab$src[i], " ", tab$call[i], "\n", sep = "")
}

# ---------------------------------------------------- hypothesis strings

bex <- brms:::rename_pars(get("brmsfit_example1", envir = asNamespace("brms")))
brms_grammar <- function(s) {
  tryCatch({
    brms::hypothesis(bex, s)
    TRUE
  }, error = function(e) {
    # brms's grammar refusal, told from its other errors by its text; a
    # refusal of an unknown PARAMETER means the grammar was accepted
    !grepl("Every hypothesis must be of the form", conditionMessage(e),
           fixed = TRUE)
  })
}
frm_grammar <- function(s) {
  tryCatch({
    frmtmb:::hyp_parse(s)
    TRUE
  }, error = function(e) FALSE)
}

hs <- unique(unlist(lapply(hyp_hits, function(x) x$h)))
hsrc <- vapply(hs, function(s) {
  i <- which(vapply(hyp_hits, function(x) s %in% x$h, NA))[1L]
  hyp_hits[[i]]$src
}, "")
hb <- vapply(hs, brms_grammar, NA)
hf <- vapply(hs, frm_grammar, NA)
cat("\n== D3 hypothesis grammar, distinct strings:", length(hs), "\n")
cat("  both accept :", sum(hb & hf), "\n")
cat("  both refuse :", sum(!hb & !hf), "\n")
cat("  FALSE ALARM (frmtmb refuses, brms accepts):", sum(hb & !hf), "\n")
for (i in which(hb & !hf)) cat("    ", hsrc[i], " '", hs[i], "'\n", sep = "")
cat("  MISSED (frmtmb accepts, brms refuses):", sum(!hb & hf), "\n")
for (i in which(!hb & hf)) cat("    ", hsrc[i], " '", hs[i], "'\n", sep = "")
cat("  refused by both, listed:\n")
for (i in which(!hb & !hf)) cat("    ", hsrc[i], " '", hs[i], "'\n", sep = "")

# The machine-readable arm, for the base-against-lane diff.
out <- file.path("dev/adefects-log", paste0("falsealarm-", tag, ".tsv"))
rows <- rbind(
  data.frame(kind = "autocor", src = tab$src, text = tab$call,
             brms = tab$brms, frm = tab$frm),
  data.frame(kind = "hypothesis", src = unname(hsrc), text = hs,
             brms = unname(hb), frm = unname(hf)))
utils::write.table(rows, out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("wrote", out, nrow(rows), "rows\n")

cat("\nDONE\n")
