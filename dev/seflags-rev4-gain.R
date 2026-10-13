# Reviewer of lane seflags, check of punch 2b: the ill warning's
# "stopped short" report. Degree-8 raw polynomial through I(x^k), with
# and without (1 | g), seeds 1 to 20. Per fit: code, whether the ill
# warning fired, the gain it states (or none), and the true shortfall,
# the poly(x, 8) fit's log-likelihood (the same likelihood, well
# conditioned) less this fit's. A gain stated above 0.005 where the
# shortfall is below it is a false alarm; a shortfall above 0.005 with
# no gain stated is a miss.
#   Rscript dev/seflags-rev4-gain.R <lane|op4> <out.tsv>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev4-common.R")
out <- commandArgs(trailingOnly = TRUE)[2]
pw <- paste0(" + I(x^", 2:8, ")", collapse = "")
cat("re\tseed\tcode\tref_code\till\tgain\tshort\n", file = out)
for (re in c(FALSE, TRUE)) for (s in 1:20) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  ex <- if (re) " + (1 | g)" else ""
  r <- cap(frm(stats::as.formula(paste0("y ~ x", pw, ex)), data = d))
  g <- suppressWarnings(frm(stats::as.formula(paste0("y ~ poly(x, 8)", ex)),
                            data = d))
  iw <- grep("ill-conditioned", r$w, value = TRUE)
  gain <- if (length(iw) && grepl("raises the log-likelihood by", iw[1])) {
    as.numeric(sub(".*raises the log-likelihood by ([0-9.e+-]+),.*", "\\1",
                   iw[1]))
  } else NA
  short <- as.numeric(logLik(g)) - as.numeric(logLik(r$v))
  cat(re, s, r$v$opt$convergence, g$opt$convergence, length(iw) > 0, gain,
      signif(short, 4), sep = "\t", file = out, append = TRUE)
  cat("\n", file = out, append = TRUE)
}
cat("done\n")
