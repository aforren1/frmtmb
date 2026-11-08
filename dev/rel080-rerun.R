# Run a lane's or a reviewer's script on the release build: every lane
# library (wt-<lane>-lib, rev-<lane>-lib) in it, and in each file it
# sources, reads rellib-r19, and every path into a lane tree reads the
# release tree, so that nothing is written into a lane tree. rellib-r18
# stays where the script puts it (behind rellib-r19, or a before arm).
#
#   Rscript dev/rel080-rerun.R <script> [args to the script ...]
#
# from the release tree's root. The script's commandArgs(TRUE) are the
# args after <script>.
local({
  a <- commandArgs(trailingOnly = TRUE)
  REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel080"
  fix <- function(s) {
    s <- gsub("C:/Users/adf44/source/r/(wt|rev)-[a-z0-9]+-lib",
              "C:/Users/adf44/source/r/rellib-r19", s)
    s <- gsub("C:/Users/adf44/source/r/frmtmb-wt-(ci|latperf|mfxse|loglik|osa|nestquad|arstep|fitedge|bigx|nlstart|readers|superassign)",
              REL, s)
    s <- gsub("C:\\\\Users\\\\adf44\\\\source\\\\r\\\\frmtmb-wt-[a-z]+",
              "C:\\\\Users\\\\adf44\\\\source\\\\r\\\\frmtmb-wt-rel080", s)
    s
  }
  rest <- a[-1L]
  script <- a[1L]
  assign("commandArgs", function(trailingOnly = FALSE) {
    if (trailingOnly) rest else c("Rscript", paste0("--file=", script), "--args", rest)
  }, envir = globalenv())
  assign("source", function(file, local = FALSE, ...) {
    env <- if (isTRUE(local)) parent.frame() else if (is.environment(local)) local else globalenv()
    invisible(eval(parse(text = fix(readLines(file, warn = FALSE)),
                         keep.source = FALSE), envir = env))
  }, envir = globalenv())
  setwd(REL)
  cat("rerun:", script, paste(rest, collapse = " "), "\n")
  eval(parse(text = fix(readLines(script, warn = FALSE)), keep.source = FALSE),
       envir = globalenv())
})
