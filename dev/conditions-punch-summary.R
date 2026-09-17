# Lane wt-conditions, punch round 1: the generated blocks of the "Punch
# round 1" section of dev/conditions-findings.md, read from a sweep made
# by dev/conditions-rev-runset.sh (one test file per R process, runner
# dev/conditions-rev-sweep-onefile.R, which saves each caught condition
# as a row of an .rds file, so no text is parsed back).
#   Rscript dev/conditions-punch-summary.R <sweep dir>
d <- commandArgs(trailingOnly = TRUE)[1L]
num_after <- function(x, k) {
  as.integer(vapply(regmatches(x, regexec(paste0(k, " +([0-9]+)"), x)),
                    function(m) if (length(m) > 1L) m[2L] else NA_character_,
                    ""))
}
s <- readLines(file.path(d, "summary.txt"), warn = FALSE)
files <- c(Sys.glob("tests/testthat/test-*.R"),
           Sys.glob("extensions/*/tests/testthat/test-*.R"))
cat("== suites, one file per R process ==\n")
cat("test files on disk", length(files), " result lines", length(s), "\n")
pkg_of <- function(x) ifelse(grepl("extensions/", x),
                             sub("^.*extensions/([^/]+)/.*$", "\\1", x),
                             "frmtmb")
pkg <- pkg_of(s)
for (p in sort(unique(pkg_of(files)))) {
  x <- s[pkg == p]
  ok <- grepl("^BLOCKS", x)
  cat(sprintf(paste0("%-16s files %3d of %3d (no result line %d) blocks %4d",
                     " PASS %5d FAIL %d ERROR %d SKIP %d\n"),
              p, sum(ok), sum(pkg_of(files) == p), sum(!ok),
              sum(num_after(x[ok], "BLOCKS")), sum(num_after(x[ok], "PASS")),
              sum(num_after(x[ok], "FAIL")), sum(num_after(x[ok], "ERROR")),
              sum(num_after(x[ok], "SKIP"))))
}
bad <- s[!grepl("^BLOCKS", s) | num_after(s, "FAIL") > 0 |
           num_after(s, "ERROR") > 0]
if (length(bad)) cat("not clean:\n", paste0("  ", bad, "\n"), sep = "")
cat("start", readLines(file.path(d, ".start")), " end",
    readLines(file.path(d, ".end")), "\n")
# which library the sampler suite loaded Stan from: a fresh compile
# against the wrong StanHeaders fails in compileCode()
logs <- list.files(d, "frmtmb[.]sample.*[.]log$", full.names = TRUE)
txt <- unlist(lapply(logs, readLines, warn = FALSE))
cat("frmtmb.sample logs", length(logs), " lines naming compileCode",
    sum(grepl("compileCode", txt)), "\n")

rows <- do.call(rbind, lapply(list.files(d, "[.]rds$", full.names = TRUE),
                              readRDS))
err <- rows[grepl("^expect_error ", rows$key), ]
err$pkg <- pkg_of(err$file)
cat("\n== expect_error() conditions ==\n")
tab <- t(vapply(split(err, err$pkg), function(y) {
  c(caught = nrow(y), frmtmb_error = sum(grepl("frmtmb_error", y$class)),
    other = sum(!grepl("frmtmb_error", y$class)))
}, numeric(3)))
print(rbind(tab, ALL = colSums(tab)))

cat("\n== extension subclass, in each extension's own suite ==\n")
cat("(the construction of dev/conditions-rev-subclass-sweep.R)\n")
ext <- err[err$pkg != "frmtmb" & grepl("frmtmb_error", err$class), ]
sub_tab <- t(vapply(split(ext, ext$pkg), function(y) {
  own <- paste0(gsub(".", "_", y$pkg[1L], fixed = TRUE), "_error")
  has <- vapply(strsplit(y$class, "/"), function(k) own %in% k, NA)
  c(caught = nrow(y), with_own = sum(has), without = sum(!has))
}, numeric(3)))
sub_tab <- rbind(sub_tab, ALL = colSums(sub_tab))
print(cbind(sub_tab, pct_with = round(100 * sub_tab[, 2] / sub_tab[, 1], 1)))
for (p in sort(unique(ext$pkg))) {
  own <- paste0(gsub(".", "_", p, fixed = TRUE), "_error")
  y <- ext[ext$pkg == p, ]
  no <- y[!vapply(strsplit(y$class, "/"), function(k) own %in% k, NA), ]
  if (!nrow(no)) next
  cat("\n", p, ": ", nrow(no), " without ", own, "\n", sep = "")
  u <- unique(data.frame(cls = sub("/frmtmb_error.*$|/error.*$", "",
                                   no$class),
                         msg = substr(gsub("[[:space:]]+", " ", no$message),
                                      1, 88)))
  for (i in seq_len(nrow(u))) cat("  [", u$cls[i], "] ", u$msg[i], "\n",
                                  sep = "")
}

cat("\n== expect_error() conditions that are not a frmtmb_error ==\n")
o <- err[!grepl("frmtmb_error", err$class), ]
for (i in seq_len(nrow(o))) {
  cat(sprintf("[%d] %s\n    %s\n    %s | %s\n", i,
              sub("^.*tests/testthat/", "", o$file[i]),
              substr(sub("^expect_error ", "", o$key[i]), 1, 110),
              sub("/error/condition$", "", o$class[i]),
              substr(gsub("\\\\n", " ", o$message[i]), 1, 90)))
}
