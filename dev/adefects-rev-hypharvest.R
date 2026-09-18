# Independent harvest of literal hypothesis() strings from a tree, and
# the verdict of brms, of the frmtmb build on the path, on each.
source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages({ library(frmtmb); library(brms) })
TREE <- Sys.getenv("REV_TREE")
OUT <- Sys.getenv("REV_OUT")
cat("tree", TREE, "\n")

lits <- new.env(parent = emptyenv())
add <- function(s, file) {
  if (!is.character(s) || length(s) != 1L || is.na(s)) return()
  k <- s
  cur <- get0(k, envir = lits)
  assign(k, unique(c(cur, file)), envir = lits)
}

walk <- function(e, file) {
  if (is.call(e)) {
    fn <- e[[1L]]
    nm <- if (is.name(fn)) as.character(fn) else
      if (is.call(fn) && length(fn) == 3L &&
          as.character(fn[[1L]]) %in% c("::", ":::"))
        as.character(fn[[3L]]) else ""
    if (nm == "hypothesis" && length(e) >= 2L) {
      a <- as.list(e)[-1L]
      an <- names(a); if (is.null(an)) an <- rep("", length(a))
      # brms's second positional argument, or the named one
      idx <- which(an == "hypothesis")
      if (!length(idx)) {
        pos <- which(an == "")
        idx <- if (length(pos) >= 2L) pos[2L] else integer(0)
      }
      if (length(idx)) {
        h <- a[[idx[1L]]]
        if (is.character(h)) for (s in h) add(s, file)
        else if (is.call(h) && identical(h[[1L]], as.name("c"))) {
          for (el in as.list(h)[-1L]) if (is.character(el)) add(el, file)
        } else if (is.call(h) && identical(h[[1L]], as.name("paste0"))) {
          # paste0(<expr>, " = 0") and friends: not a literal, skipped
        }
      }
    }
    for (i in seq_along(e)) if (!is.null(e[[i]])) walk(e[[i]], file)
  } else if (is.pairlist(e) || is.expression(e) || is.list(e)) {
    for (i in seq_along(e)) if (!is.null(e[[i]])) walk(e[[i]], file)
  }
}

files <- c(
  list.files(file.path(TREE, "R"), "[.]R$", full.names = TRUE),
  list.files(file.path(TREE, "tests"), "[.]R$", full.names = TRUE,
             recursive = TRUE),
  list.files(file.path(TREE, "vignettes"), "[.]Rmd$", full.names = TRUE,
             recursive = TRUE),
  list.files(file.path(TREE, "extensions"), "[.](R|Rmd)$",
             full.names = TRUE, recursive = TRUE),
  list.files(file.path(TREE, "man"), "[.]Rd$", full.names = TRUE))
files <- files[!grepl("test-adefects[.]R$", files)]

for (f in files) {
  txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) NULL)
  if (is.null(txt)) next
  if (grepl("[.]Rmd$", f)) {
    keep <- rep(FALSE, length(txt)); inch <- FALSE
    for (i in seq_along(txt)) {
      if (grepl("^```", txt[i])) { inch <- !inch; next }
      keep[i] <- inch
    }
    txt <- txt[keep]
  }
  if (grepl("[.]Rd$", f)) {
    # examples and usage blocks only; parse what parses
    txt <- txt[!grepl("^\\\\", txt)]
  }
  ex <- tryCatch(parse(text = txt, keep.source = FALSE),
                 error = function(e) NULL)
  if (is.null(ex)) next
  walk(ex, basename(f))
}

hs <- ls(lits)
cat("distinct literal hypothesis strings:", length(hs), "\n")

bfit <- get("brmsfit_example1", envir = asNamespace("brms"))
brms_verdict <- function(h) {
  r <- tryCatch({ brms::hypothesis(bfit, h); "accept" },
                error = function(e) {
                  if (grepl("Every hypothesis must be of the form",
                            conditionMessage(e), fixed = TRUE))
                    "grammar-refuse" else "accept"
                })
  r
}
frm_verdict <- function(h) {
  tryCatch({ frmtmb:::hyp_parse(h); "accept" },
           error = function(e) {
             m <- conditionMessage(e)
             if (grepl("states no relation", m, fixed = TRUE))
               "d3-refuse" else "other-refuse"
           })
}

res <- data.frame(h = hs,
                  brms = vapply(hs, brms_verdict, ""),
                  frm = vapply(hs, frm_verdict, ""),
                  files = vapply(hs, function(k)
                    paste(get(k, envir = lits), collapse = ";"), ""),
                  stringsAsFactors = FALSE)
write.table(res, OUT, sep = "\t", row.names = FALSE, quote = FALSE)
cat("brms accept:", sum(res$brms == "accept"),
    " brms grammar-refuse:", sum(res$brms == "grammar-refuse"), "\n")
cat("frm accept:", sum(res$frm == "accept"),
    " frm d3-refuse:", sum(res$frm == "d3-refuse"),
    " frm other-refuse:", sum(res$frm == "other-refuse"), "\n")
cat("FALSE ALARM (frm refuses, brms accepts):",
    sum(res$frm != "accept" & res$brms == "accept"), "\n")
cat("MISS (frm accepts, brms grammar-refuses):",
    sum(res$frm == "accept" & res$brms == "grammar-refuse"), "\n")
bad <- res[res$frm != "accept" & res$brms == "accept", ]
if (nrow(bad)) print(bad[, c("h", "frm", "files")])
mis <- res[res$frm == "accept" & res$brms == "grammar-refuse", ]
if (nrow(mis)) print(mis[, c("h", "files")])
