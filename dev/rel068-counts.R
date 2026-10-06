# Every count dev/round-20261005.md quotes, generated from the logs of
# the 0.68.0 release tree, so none is typed.
#
#   Rscript dev/rel068-counts.R > dev/rel068-log/counts.md
#
# The suite is every test file of the eight packages with every gate
# set (dev/release/suite-files/, one R process per file). The files the
# ledger regeneration rewrote, and test-gp-by.R after the runner fix,
# were run again after it (dev/release/tier-files/); their later result
# replaces the earlier one, and both are listed.
rel <- "dev/release"
jobs <- read.table("dev/rel068-log/suite.jobs", col.names = c("pkg", "path"),
                   stringsAsFactors = FALSE)
parse_res <- function(f) {
  if (!file.exists(f)) return(NULL)
  x <- readLines(f, warn = FALSE)
  r <- grep("^RESULT ", x, value = TRUE)
  if (!length(r)) return(NULL)
  r <- r[length(r)]
  kv <- regmatches(r, gregexpr("[a-z]+=[0-9]+", r))[[1]]
  v <- as.integer(sub("^[a-z]+=", "", kv))
  names(v) <- sub("=.*$", "", kv)
  list(v = v, lib = grep("^lib: ", x, value = TRUE)[1],
       skips = sub("^DETAIL SKIP ", "", grep("^DETAIL SKIP ", x, value = TRUE)))
}
rows <- list()
reruns <- character()
for (i in seq_len(nrow(jobs))) {
  b <- paste0(jobs$pkg[i], "--", sub("[.]R$", "", basename(jobs$path[i])),
              ".txt")
  first <- parse_res(file.path(rel, "suite-files", b))
  later <- parse_res(file.path(rel, "tier-files", b))
  use <- if (!is.null(later)) later else first
  if (!is.null(later)) {
    reruns <- c(reruns, sprintf("%s `%s`: first pass %d fail %d err %d, rerun pass %d fail %d err %d",
                                jobs$pkg[i], basename(jobs$path[i]),
                                first$v[["pass"]], first$v[["fail"]],
                                first$v[["err"]], later$v[["pass"]],
                                later$v[["fail"]], later$v[["err"]]))
  }
  rows[[i]] <- data.frame(
    pkg = jobs$pkg[i], file = basename(jobs$path[i]),
    ran = !is.null(use),
    pass = if (is.null(use)) NA else use$v[["pass"]],
    fail = if (is.null(use)) NA else use$v[["fail"]],
    err = if (is.null(use)) NA else use$v[["err"]],
    skip = if (is.null(use)) NA else use$v[["skip"]],
    warn = if (is.null(use)) NA else use$v[["warn"]],
    rellib = !is.null(use) && grepl("rellib-r6", use$lib) &&
      !grepl("rellib-r5|wt-", use$lib),
    skips = if (is.null(use)) "" else paste(use$skips, collapse = " | "))
}
d <- do.call(rbind, rows)
cat("<!-- BEGIN GENERATED: dev/rel068-counts.R -->\n")
cat("\n### All eight suites, every gate set (brms, drmTMB, fuzz)\n\n")
cat("From `dev/release/suite-files/` (", sum(d$ran), " of ", nrow(d),
    " files with a RESULT line; ", sum(d$rellib), " whose `lib:` line ",
    "names rellib-r6 for the package and for frmtmb).\n\n", sep = "")
cat("| package | files | pass | fail | error | skip | warn |\n")
cat("|---|---|---|---|---|---|---|\n")
for (p in unique(d$pkg)) {
  s <- d[d$pkg == p, ]
  cat(sprintf("| %s | %d | %d | %d | %d | %d | %d |\n", p, nrow(s),
              sum(s$pass), sum(s$fail), sum(s$err), sum(s$skip),
              sum(s$warn)))
}
tot <- c(nrow(d), sum(d$pass), sum(d$fail), sum(d$err), sum(d$skip),
         sum(d$warn))
cat("| ", paste0("**", c("total", tot), "**", collapse = " | "), " |\n",
    sep = "")
sk <- d[d$skip > 0, ]
cat("\nEvery skip, by file and test:\n\n")
for (i in seq_len(nrow(sk))) {
  cat("- ", sk$pkg[i], " `", sk$file[i], "` (", sk$skip[i], "): ",
      sk$skips[i], "\n", sep = "")
}
cat("\nFiles run again after the ledger regeneration and the runner fix ",
    "(`dev/release/tier-files/`):\n\n", sep = "")
for (r in reruns) cat("- ", r, "\n", sep = "")
bad <- d[d$fail > 0 | d$err > 0 | d$warn > 0, ]
cat("\nFiles with a failure, an error or an escaped warning after the ",
    "reruns: ", nrow(bad), if (nrow(bad)) paste0(" (",
    paste(bad$file, collapse = ", "), ")"), ".\n", sep = "")

cat("\n### Scale tier (FRMTMB_SCALE_TESTS set)\n\nFrom `dev/release/scale.log`. ")
sl <- readLines(file.path(rel, "scale.log"), warn = FALSE)
sr <- grep("^RESULT ", sl, value = TRUE)
snum <- function(k) sum(as.integer(sub(paste0(".* ", k, "=([0-9]+).*"),
                                       "\\1", sr)))
cat(grep(" ran [0-9]+ of [0-9]+$", sl, value = TRUE), "; pass ", snum("pass"),
    ", fail ", snum("fail"), ", error ", snum("err"), ", skip ", snum("skip"),
    ", warn ", snum("warn"), ".\n", sep = "")

cat("\n### The ledger (`dev/brmsport-ledger.tsv`)\n\n")
lb <- utils::read.delim("dev/rel068-ledger-before.tsv", quote = "",
                        colClasses = "character")
la <- utils::read.delim("dev/brmsport-ledger.tsv", quote = "",
                        colClasses = "character")
lev <- unique(c(lb$outcome, la$outcome))
ob <- table(factor(lb$outcome, levels = lev))
oa <- table(factor(la$outcome, levels = lev))
cat("| outcome | 0.67.0 | 0.68.0 |\n|---|---|---|\n")
for (k in lev) cat(sprintf("| %s | %d | %d |\n", k, ob[[k]], oa[[k]]))
cat(sprintf("| **total** | **%d** | **%d** |\n", nrow(lb), nrow(la)))
key <- function(x) paste(x$file, x$line)
m <- match(key(la), key(lb))
moved <- which(la$outcome != lb$outcome[m] | la$class != lb$class[m] |
                 la$reason != lb$reason[m])
cat("\nBin 1 passes: ", sum(lb$outcome == "pass"), " of ", nrow(lb),
    " before, ", sum(la$outcome == "pass"), " of ", nrow(la),
    " after.\nRows whose outcome, class or reason moved: ", length(moved),
    "\n\n", sep = "")
for (i in moved) {
  cat(sprintf("- `%s`: %s / %s -> %s / %s\n", la$id[i], lb$outcome[m[i]],
              if (nzchar(lb$class[m[i]])) lb$class[m[i]] else "-",
              la$outcome[i], if (nzchar(la$class[i])) la$class[i] else "-"))
}

cat("\n### Ported brms suite, verdicts asserted\n\n",
    "From `dev/rel068-log/tier.txt`. ", sep = "")
tl <- readLines("dev/rel068-log/tier.txt", warn = FALSE)
tr <- grep("^RESULT test-brms-suite", tl, value = TRUE)
num <- function(k) sum(as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1",
                                      tr)))
cat(length(tr), " files; pass ", num("pass"), ", fail ", num("fail"),
    ", error ", num("err"), ", skip ", num("skip"), ".\n", sep = "")

ws <- file.path(rel, "warnscan-068")
cat("\n### Escaped warnings (`dev/release/warnscan-068/`)\n\n")
if (file.exists(file.path(ws, "scanned.tsv"))) {
  sc <- utils::read.delim(file.path(ws, "scanned.tsv"), header = FALSE,
                          col.names = c("tier", "file", "n"))
  nof <- if (file.exists(file.path(ws, "nofile.txt"))) {
    readLines(file.path(ws, "nofile.txt"), warn = FALSE)
  } else character()
  esc <- readLines(file.path(ws, "escaped.tsv"), warn = FALSE)
  cat("- every gate set: ", nrow(sc), " files scanned; ",
      sum(sc$n > 0) - length(nof), " with an escaped warning; ",
      length(nof), " the scan could not run\n", sep = "")
  cat("- escaped warnings: ", length(esc), "\n", sep = "")
} else {
  cat("No scan.\n")
}

cat("\n### R CMD check --as-cran\n\n",
    "From `dev/release/check-068/<pkg>/check.log`.\n\n", sep = "")
for (dd in list.dirs(file.path(rel, "check-068"), recursive = FALSE)) {
  f <- file.path(dd, "check.log")
  st <- if (file.exists(f)) grep("^Status:", readLines(f, warn = FALSE),
                                 value = TRUE) else "no log"
  nt <- if (file.exists(f)) {
    x <- readLines(f, warn = FALSE)
    x[grepl("\\.\\.\\. (\\[[^]]*\\] )?(NOTE|WARNING|ERROR)$", x)]
  } else character()
  cat("- ", basename(dd), ": ", if (length(st)) st else "no Status line",
      if (length(nt)) paste0(" (", paste(trimws(nt), collapse = "; "), ")"),
      "\n", sep = "")
}
cat("<!-- END GENERATED -->\n")
