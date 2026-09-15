root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (length(gregexpr(old, txt, fixed = TRUE)[[1]]) != 1L)
    stop("many in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

# The reviewer's library still holds the SWAP build of this same lane
# (dev/generics-whichbuild.R confirms which mechanism each library
# carries), so the two designs can be timed in ONE interleaved run
# instead of compared across runs on different machine load.
sub1("dev/generics-loadcost.R",
paste0("FIX <- \"C:/Users/adf44/source/r/generics-lib\"\n",
       "BASE <- \"C:/Users/adf44/source/r/rellib-r3\"\n"),
paste0("FIX <- \"C:/Users/adf44/source/r/generics-lib\"\n",
       "BASE <- \"C:/Users/adf44/source/r/rellib-r3\"\n",
       "# the review's library, which holds the SWAP build of this lane.\n",
       "# READ ONLY; nothing is ever installed there.\n",
       "SWAP <- \"C:/Users/adf44/source/r/genrev-lib\"\n"))

sub1("dev/generics-loadcost.R",
paste0("  `frmtmb FIX`  = list(FIX,  \"library(frmtmb)\"),\n",
       "  `frmtmb BASE` = list(BASE, \"library(frmtmb)\"),\n"),
paste0("  `frmtmb FIX (active)` = list(FIX,  \"library(frmtmb)\"),\n",
       "  `frmtmb SWAP (round 1)` = list(SWAP, \"library(frmtmb)\"),\n",
       "  `frmtmb BASE` = list(BASE, \"library(frmtmb)\"),\n"))

sub1("dev/generics-loadcost.R",
'pp <- perm_p(hr[, "frmtmb FIX"], hr[, "frmtmb BASE"], nperm)',
paste0('pp <- perm_p(hr[, "frmtmb FIX (active)"], hr[, "frmtmb BASE"],\n',
       '             nperm)\n',
       'sw <- perm_p(hr[, "frmtmb SWAP (round 1)"], hr[, "frmtmb BASE"],\n',
       '             nperm)\n',
       'ab <- perm_p(hr[, "frmtmb FIX (active)"],\n',
       '             hr[, "frmtmb SWAP (round 1)"], nperm)'))

sub1("dev/generics-loadcost.R",
paste0('cat(sprintf("\\nFIX - BASE on minima, hi-res:     %+.4f s  p = %.4f\\n",\n',
       '            pp[["obs"]], pp[["p"]]))\n'),
paste0('cat(sprintf("\\nACTIVE - BASE on minima, hi-res:  %+.4f s  p = %.4f\\n",\n',
       '            pp[["obs"]], pp[["p"]]))\n',
       'cat(sprintf("SWAP   - BASE, same:              %+.4f s  p = %.4f\\n",\n',
       '            sw[["obs"]], sw[["p"]]))\n',
       'cat(sprintf("ACTIVE - SWAP, same:              %+.4f s  p = %.4f\\n",\n',
       '            ab[["obs"]], ab[["p"]]))\n'))

sub1("dev/generics-loadcost.R",
paste0('cat(sprintf("FIX - BASE on proc.time() minima: %+.4f s (%.1f ticks)\\n",\n',
       '            min(st[, "frmtmb FIX"]) - min(st[, "frmtmb BASE"]),\n',
       '            (min(st[, "frmtmb FIX"]) - min(st[, "frmtmb BASE"])) /\n',
       '              pt[["min"]]))\n'),
paste0('cat(sprintf("ACTIVE - BASE on proc.time():      %+.4f s (%.1f ticks)\\n",\n',
       '            min(st[, "frmtmb FIX (active)"]) -\n',
       '              min(st[, "frmtmb BASE"]),\n',
       '            (min(st[, "frmtmb FIX (active)"]) -\n',
       '               min(st[, "frmtmb BASE"])) / pt[["min"]]))\n'))

sub1("dev/generics-loadcost.R",
paste0('cat(sprintf("  FIX minus that arm:                 %+.4f s\\n",\n',
       '            min(hr[, "frmtmb FIX"]) - min(hr[, "BASE + nlme + gen"])))\n',
       'cat("  so the two new hard Imports are the whole of it.\\n")\n'),
paste0('cat(sprintf("  ACTIVE minus that arm:              %+.4f s\\n",\n',
       '            min(hr[, "frmtmb FIX (active)"]) -\n',
       '              min(hr[, "BASE + nlme + gen"])))\n',
       'cat(sprintf("  SWAP   minus that arm:              %+.4f s\\n",\n',
       '            min(hr[, "frmtmb SWAP (round 1)"]) -\n',
       '              min(hr[, "BASE + nlme + gen"])))\n'))
cat("DONE\n")
