# Lane eamhier: where the 2.4 GB of a 30 x 400 hierarchical Wiener fit
# actually goes.
#
# The replicate processes peak near 2.4 GB of process working set while
# R's own heap stays at 266 MB, so the difference is outside R and the
# question is whether it is the DATA (12,000 rows is nothing) or the
# RTMB tape. This walks the stages and reads the process working set
# from the operating system at each one, because gc() cannot see memory
# R did not allocate.
#
# The reading is the WORKING SET, which the operating system may trim
# under pressure, so each stage reports the peak working set as well.
#
# Run: Rscript --vanilla dev/eamhier-scripts/eamhier-memory.R
source("dev/eamhier-scripts/eamhier-common.R")
eamhier_libs()
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

ws <- function() {
  cmd <- sprintf(
    "$p=Get-Process -Id %d; '{0} {1}' -f $p.WorkingSet64,$p.PeakWorkingSet64",
    Sys.getpid())
  out <- suppressWarnings(system2(
    "powershell", c("-NoProfile", "-Command", cmd), stdout = TRUE))
  as.numeric(strsplit(trimws(out[length(out)]), " ")[[1L]]) / 1024^2
}
# to stderr, which R does not buffer: the first run of this script was
# redirected to a file and read back with its first six stages missing.
stage <- function(tag) {
  m <- ws()
  cat(sprintf("%-34s ws %8.1f MB  peak %8.1f MB  R heap %7.1f MB\n",
              tag, m[1L], m[2L], sum(gc()[, 2L])), file = stderr())
  invisible(m)
}

# "prefit" stops before the fit, which is the expensive half in TIME and
# does not move the peak.
only_prefit <- identical(commandArgs(trailingOnly = TRUE)[1L],
                         "prefit")

stage("R plus the two packages")
d <- eamhier_data(20260908L, "A")
stage(paste0("the data, ", nrow(d), " rows"))
form <- eamhier_form("A", "pg")
fr <- frm(form, family = wiener(), data = d, dry_run = "frame")
stage("frm(dry_run = 'frame')")
ob <- frm(form, family = wiener(), data = d, dry_run = "objective")
stage("frm(dry_run = 'objective'), the tape")
v <- ob$obj$fn(ob$obj$par)
stage("one objective evaluation")
g <- ob$obj$gr(ob$obj$par)
stage("one gradient")
for (i in 1:5) g <- ob$obj$gr(ob$obj$par)
stage("five more gradients")
rm(ob, fr)
invisible(gc())
stage("after dropping the tape and gc()")
if (only_prefit) quit(save = "no")
fit <- frm(form, family = wiener(), data = d, se = FALSE)
stage("the whole fit, se = FALSE")
# confint() is the public route that triggers sdreport(), which is
# where a fit with se = TRUE pays for its standard errors.
ci <- suppressWarnings(stats::confint(fit))
stage("confint(), which runs sdreport()")
