# Item 1562: which count each information-criterion path reads on a
# multivariate fit whose subset() leaves rows to no response.
# Usage: Rscript dev/silent-nobs-probe.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
}
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")

# the triage's construction (round-briefs/triage-081/scratch-2, run H)
set.seed(5)
d6 <- data.frame(x = rnorm(80))
d6$y1 <- d6$x + rnorm(80)
d6$y2 <- d6$x + rnorm(80)
d6$s1 <- seq_len(80) <= 72
d6$s2 <- seq_len(80) <= 72
d6$s3 <- seq_len(80) <= 60
f <- frm(bf(y1 | subset(s1) ~ x) + bf(y2 | subset(s2) ~ x) +
           set_rescor(FALSE), data = d6)
ll <- logLik(f)
cat("A: both responses rows 1..72 of 80\n")
cat("  logLik nobs attr", attr(ll, "nobs"), " nobs()", nobs(f),
    " nobs(resp=y1)", nobs(f, resp = "y1"), "\n")
cat("  BIC", format(BIC(f), digits = 12),
    " -2ll + log(72) df", format(-2 * as.numeric(ll) +
                                   log(72) * attr(ll, "df"), digits = 12),
    " -2ll + log(80) df", format(-2 * as.numeric(ll) +
                                   log(80) * attr(ll, "df"), digits = 12),
    "\n")
cat("  df.residual", df.residual(f), " df", attr(ll, "df"), "\n")
lk <- tryCatch(log_lik(f), error = function(e) conditionMessage(e))
cat("  log_lik dim", if (is.matrix(lk)) dim(lk) else lk, "\n")
lo <- tryCatch(loo(f), error = function(e) conditionMessage(e))
cat("  loo pointwise rows",
    if (is.list(lo)) nrow(lo$pointwise) else lo, "\n")

# different subsets: y1 rows 1..72, y2 rows 1..60; union 72
g <- frm(bf(y1 | subset(s1) ~ x) + bf(y2 | subset(s3) ~ x) +
           set_rescor(FALSE), data = d6)
cat("B: y1 rows 1..72, y2 rows 1..60\n")
cat("  logLik nobs attr", attr(logLik(g), "nobs"), " nobs()", nobs(g),
    " df.residual", df.residual(g), "\n")

# the same model on all 80 rows: anova of it against A compares
# likelihoods of different data
h <- frm(bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(FALSE), data = d6)
h2 <- frm(bf(y1 | subset(s1) ~ 1) + bf(y2 | subset(s2) ~ x) +
            set_rescor(FALSE), data = d6)
cat("C: anova(80-row fit, 72-row fit):",
    tryCatch({a <- anova(h2, h); "ran"},
             error = function(e) conditionMessage(e)), "\n")
cat("  anova(72-row nested fits):",
    tryCatch({a <- anova(h2, f); paste("ran, Chisq", format(a$Chisq[2]))},
             error = function(e) conditionMessage(e)), "\n")
# the absent case: no subset, all 80 rows
cat("D: no subset: logLik nobs attr", attr(logLik(h), "nobs"),
    " nobs()", nobs(h), " df.residual", df.residual(h), "\n")
