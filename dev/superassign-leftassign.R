# pkgcheck's left-assign check, run alone.
#
# pkgcheck 0.1.x decides "Package uses global assignment operator
# ('<<-')" in pkgchk_left_assign() (R/checks-left-assign.R upstream):
# every *.R/*.r/*.q/*.s file under R/ is parsed, and the LEFT_ASSIGN and
# EQ_ASSIGN tokens of utils::getParseData() are tallied by text. The
# check fails when the `<<-` tally is above zero after two exemptions
# (a `<<-` inside a setRefClass() definition, and `<<- memoise` inside
# .onLoad), neither of which frmtmb or an extension uses. The same
# function also scores "inconsistent assignment operators": it fails
# when the package uses both `<-` and `=` as an assignment.
#
# This script repeats that tally without pkgcheck, which is not
# installed here and takes hours to run in full. It also counts `->>`
# (RIGHT_ASSIGN), which pkgcheck does not, because it is the same
# operator written the other way round.
#
# Usage: Rscript dev/superassign-leftassign.R <pkgdir> [<pkgdir> ...]

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- "."

tally_dir <- function(pkg) {
  rdir <- file.path(pkg, "R")
  fl <- list.files(rdir, pattern = "\\.(q|r|s)$", ignore.case = TRUE,
                   full.names = TRUE)
  res <- lapply(fl, function(f) {
    p <- utils::getParseData(parse(f, keep.source = TRUE))
    la <- p$text[p$token %in% c("EQ_ASSIGN", "LEFT_ASSIGN")]
    ra <- p$text[p$token == "RIGHT_ASSIGN"]
    data.frame(file = basename(f),
               left = sum(la == "<-"), eq = sum(la == "="),
               global = sum(la == "<<-"), rglobal = sum(ra == "->>"),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

for (pkg in args) {
  tab <- tally_dir(pkg)
  hit <- tab[tab$global > 0L | tab$rglobal > 0L, , drop = FALSE]
  cat(sprintf("PACKAGE %s files=%d\n", basename(normalizePath(pkg)),
              nrow(tab)))
  for (i in seq_len(nrow(hit))) {
    cat(sprintf("  %s: <<- %d, ->> %d\n", hit$file[i], hit$global[i],
                hit$rglobal[i]))
  }
  cat(sprintf(paste0("  TOTAL <<- %d in %d files; ->> %d; ",
                     "'<-' %d; '=' %d\n"),
              sum(tab$global), sum(tab$global > 0L), sum(tab$rglobal),
              sum(tab$left), sum(tab$eq)))
  cat(sprintf("  pkgcheck global_assign: %s\n",
              if (sum(tab$global) > 0L) "FAIL" else "pass"))
  cat(sprintf("  pkgcheck left_assign consistency: %s\n",
              if (sum(tab$left) == 0L || sum(tab$eq) == 0L) "pass" else
                "FAIL"))
}
