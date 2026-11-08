# Item 3: how far the marginal route's standard errors are from the z tape's,
# against a bound on the condition number of the marginal covariance,
# kappa <= 1 + u sd^2 / sigma^2 (u distinct positions; a correlation's
# diagonal is 1, so its largest eigenvalue is at most u, and the smallest of V
# is at least sigma^2). Reads the runs of dev/fitedge-se-verdict.R with the
# tag given (force: FE_FORCE=1; lane: the shipped gate).
o <- "C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-out"
tag <- commandArgs(TRUE)[1]
fs <- list.files(o, sprintf("^severdict-.*-%s[.]tsv$", tag), full.names = TRUE)
d <- do.call(rbind, lapply(fs, read.delim, stringsAsFactors = FALSE,
                           colClasses = c(nan_old = "character",
                                          nan_new = "character")))
d$nan_old[is.na(d$nan_old)] <- ""
d$nan_new[is.na(d$nan_new)] <- ""
d$n <- suppressWarnings(as.numeric(sub(".*_n([0-9]+)_.*", "\\1", d$row)))
d$n[is.na(d$n)] <- ifelse(grepl("^r2se", d$row[is.na(d$n)]), 50, 60)
d$kappa <- 1 + d$n / d$sig_sd^2
br <- c(1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e12, Inf)
d$bin <- cut(d$kappa, br)
for (b in levels(d$bin)) {
  s <- d[d$bin %in% b, ]
  if (!nrow(s)) next
  # the weights row's marginal reads another constant, so its value is not
  # the tape's at any point
  cat(sprintf("KAPPA %-14s rows %3d route %3d valdiff max %.1e SE max rel %.2e fitted SE max rel %.2e NaN sets differ %d warnings differ %d\n",
              b, nrow(s), sum(s$route), max(s$valdiff[!grepl("weights", s$row)]),
              max(s$se_maxrel, na.rm = TRUE),
              max(s$fit_se_maxrel, na.rm = TRUE), sum(s$nan_old != s$nan_new),
              sum(!s$warn_prefix_same)))
}
cat(sprintf("TOTAL rows %d, route taken %d, NaN sets differ %d, warnings differ %d, logLik differ %d, SE max rel %.2e, fitted SE max rel %.2e\n",
            nrow(d), sum(d$route), sum(d$nan_old != d$nan_new),
            sum(!d$warn_prefix_same), sum(abs(d$ll_old - d$ll_new) > 0),
            max(d$se_maxrel, na.rm = TRUE), max(d$fit_se_maxrel, na.rm = TRUE)))
cat(sprintf("TIME z tape %.1f s, marginal %.1f s (whole frm() fits, %d rows)\n",
            sum(d$t_old), sum(d$t_new), nrow(d)))
