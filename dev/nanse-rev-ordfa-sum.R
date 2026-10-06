# Reviewer: summarize dev/nanse-rev-ordfa-log/<arm>-<model>.txt: per
# model and arm, replicates, fits with any frm() warning, with the SE
# warning, with the SE warning AND another fit-time warning, and fits
# whose SEs are not all finite.
#   Rscript dev/nanse-rev-ordfa-sum.R
d <- "C:/Users/adf44/source/r/frmtmb-wt-nanse/dev/nanse-rev-ordfa-log"
out <- NULL
for (f in list.files(d, "[.]txt$", full.names = TRUE)) {
  L <- readLines(f, warn = FALSE)
  arm <- sub("-.*", "", basename(f))
  model <- sub("^[a-z]+-", "", sub("[.]txt$", "", basename(f)))
  rep_i <- cumsum(grepl("^REP ", L))
  reps <- split(L, rep_i)[-1]
  n <- length(reps)
  fitw <- vapply(reps, function(r) {
    any(grepl("^   WARN ", r) & !grepl("^   WARN \\[se\\]", r))
  }, TRUE)
  sew <- vapply(reps, function(r) {
    any(grepl("^   WARN Standard errors are not available", r))
  }, TRUE)
  other <- vapply(reps, function(r) {
    w <- r[grepl("^   WARN ", r) & !grepl("^   WARN \\[se\\]", r)]
    any(!grepl("^   WARN Standard errors are not available", w))
  }, TRUE)
  nonfin <- vapply(reps, function(r) grepl("se_all_finite=FALSE", r[1]),
                   TRUE)
  err <- vapply(reps, function(r) grepl(" ERROR ", r[1]), TRUE)
  out <- rbind(out, data.frame(model, arm, reps = n, err = sum(err),
                               any_fit_warn = sum(fitw), se_warn = sum(sew),
                               se_and_other = sum(sew & other),
                               se_nonfinite = sum(nonfin),
                               nonfinite_silent = sum(nonfin & !fitw)))
}
out <- out[order(out$model, out$arm), ]
print(out, row.names = FALSE)
cat("\ntotals by arm:\n")
print(aggregate(cbind(reps, err, any_fit_warn, se_warn, se_and_other,
                      se_nonfinite, nonfinite_silent) ~ arm, out, sum))
