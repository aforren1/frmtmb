# Every count dev/round-20260929b.md quotes, generated from the logs of
# the 0.66.0 release tree, so none is typed.
#
#   Rscript dev/rel066-counts.R > dev/rel066-log/counts.md
rel <- "dev/release"

tier_rows <- function(log) {
  if (!file.exists(log)) return(NULL)
  x <- readLines(log, warn = FALSE)
  pkg <- NA_character_
  out <- list()
  for (l in x) {
    m <- regmatches(l, regexec("^== (\\S+) [0-9]+ files ==", l))[[1]]
    if (length(m)) { pkg <- m[2]; next }
    if (!startsWith(l, "RESULT ")) next
    f <- strsplit(l, " ", fixed = TRUE)[[1]]
    kv <- f[grepl("=", f)]
    v <- as.integer(sub("^[^=]*=", "", kv))
    names(v) <- sub("=.*$", "", kv)
    out[[length(out) + 1L]] <- data.frame(
      pkg = pkg, file = f[2], pass = v[["pass"]], fail = v[["fail"]],
      err = v[["err"]], skip = v[["skip"]],
      warn = if ("warn" %in% names(v)) v[["warn"]] else NA_integer_)
  }
  do.call(rbind, out)
}
ran_line <- function(log) grep(" ran [0-9]+ of [0-9]+$",
                               readLines(log, warn = FALSE), value = TRUE)

tier_table <- function(name, log) {
  d <- tier_rows(log)
  cat("\n### ", name, "\n\nFrom `", log, "`. ", sep = "")
  if (is.null(d)) { cat("No log.\n"); return(invisible()) }
  cat(ran_line(log), "\n\n", sep = "")
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
  if (nrow(sk)) {
    txt <- paste0("Files with a skip: ",
                  paste0(sk$pkg, " `", sk$file, "` (", sk$skip, ")",
                         collapse = ", "), ".")
    cat("\n", paste(strwrap(txt, width = 72), collapse = "\n"), "\n",
        sep = "")
  }
}

cat("<!-- BEGIN GENERATED: dev/rel066-counts.R -->\n")
tier_table("Ungated suite, all eight packages", file.path(rel, "suite.log"))
tier_table("Gated tier, with the brms, drmTMB and fuzz gates set",
           file.path(rel, "gated.log"))
tier_table("Scale tier (FRMTMB_SCALE_TESTS set)", file.path(rel, "scale.log"))
tier_table("frmtmb.ode's scale file, rerun after the wrap",
           file.path(rel, "scale-ode-rerun.log"))
tier_table("frmtmb.ode's ungated suite, rerun after the new helper",
           file.path(rel, "ode-recheck.log"))

cat("\n### Ported brms suite, verdicts asserted\n\n",
    "From `dev/rel066-log/tier.txt`. ", sep = "")
tl <- readLines("dev/rel066-log/tier.txt", warn = FALSE)
tr <- grep("^RESULT", tl, value = TRUE)
num <- function(k) sum(as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1",
                                      tr)))
cat(grep("^TIER ran", tl, value = TRUE), ";\npass ", num("pass"), ", fail ",
    num("fail"), ", error ", num("err"), ", skip ", num("skip"), ".\n",
    sep = "")

cat("\n### The ledger (`dev/brmsport-ledger.tsv`)\n\n")
lb <- utils::read.delim("dev/rel066-ledger-before.tsv", quote = "",
                        colClasses = "character")
la <- utils::read.delim("dev/brmsport-ledger.tsv", quote = "",
                        colClasses = "character")
ob <- table(factor(lb$outcome, levels = unique(c(lb$outcome, la$outcome))))
oa <- table(factor(la$outcome, levels = names(ob)))
cat("| outcome | lane defects alone | release |\n|---|---|---|\n")
for (k in names(ob)) cat(sprintf("| %s | %d | %d |\n", k, ob[[k]], oa[[k]]))
cat(sprintf("| **total** | **%d** | **%d** |\n", nrow(lb), nrow(la)))
key <- function(d) paste(d$file, d$line)
m <- match(key(la), key(lb))
moved <- which(la$outcome != lb$outcome[m] | la$class != lb$class[m])
cat("\nBin 1 passes: ", sum(lb$outcome == "pass"), " of ", nrow(lb),
    " before, ", sum(la$outcome == "pass"), " of ", nrow(la),
    " after.\nRows whose outcome or class moved: ", length(moved), ".\n",
    sep = "")

cat("\n### Escaped warnings (`dev/release/warnscan-066/`)\n\n")
esc <- file.path(rel, "warnscan-066", "escaped.tsv")
# scanned.tsv: tier, package--file, lines the scan wrote for it; a
# file the scan could not run holds the one line NOFILE, so 0 lines
# means it ran and nothing escaped
sc <- utils::read.delim(file.path(rel, "warnscan-066", "scanned.tsv"),
                        header = FALSE, col.names = c("tier", "file", "n"))
for (t in c("ungated", "gated")) {
  s <- sc[sc$tier == t, ]
  cat("- ", t, ": ", nrow(s), " files scanned, ", sum(s$n > 0),
      " with any line\n", sep = "")
}
n_esc <- if (file.exists(esc)) length(readLines(esc, warn = FALSE)) else NA
cat("- escaped warnings: ", n_esc, "\n", sep = "")

cat("\n### R CMD check --as-cran\n\n",
    "From `dev/release/check-066/<pkg>/check.log`.\n\n", sep = "")
for (d in list.dirs(file.path(rel, "check-066"), recursive = FALSE)) {
  f <- file.path(d, "check.log")
  st <- if (file.exists(f)) grep("^Status:", readLines(f, warn = FALSE),
                                 value = TRUE) else "no log"
  nt <- if (file.exists(f)) {
    x <- readLines(f, warn = FALSE)
    # a timed check line reads "... [19s] NOTE"
    x[grepl("\\.\\.\\. (\\[[^]]*\\] )?(NOTE|WARNING|ERROR)$", x)]
  } else character()
  cat("- ", basename(d), ": ", if (length(st)) st else "no Status line",
      if (length(nt)) paste0(" (", paste(trimws(nt), collapse = "; "), ")"),
      "\n", sep = "")
}
cat("<!-- END GENERATED -->\n")
