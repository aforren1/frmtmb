## What the REJECTED candidate would do at other seeds.
##
## This measures the unpaired statistic, the mean of four treatment
## ratios against the largest of four null ratios, by resampling blocks
## of four out of the 60 seeds an arm that dev/coh-calib3.R measured.
## Its answer, 4 in 20,000, is a failure rate on data that HAS the
## effect, and it is kept because it is what made the candidate look
## safe. What killed it is dev/coh-absent.R, where the same statistic
## PASSES on data that does not have the effect. The shipped assertion
## is the paired quotient, measured in dev/coh-paired.R.
##
## Run: Rscript dev/coh-falsealarm.R
kv <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[grepl("ratio=", ln, fixed = TRUE)]
  arm <- sub("^arm=([a-z]+).*$", "\\1", ln)
  ratio <- as.numeric(sub("^.*\tratio=([^\t]*).*$", "\\1", ln))
  data.frame(arm = arm, ratio = ratio)
}
d <- kv("dev/coh-calib3.tsv")
tr <- d$ratio[d$arm == "treat"]
nu <- d$ratio[d$arm == "null"]
cat("treat n=", length(tr), " min ", min(tr), " mean ", mean(tr),
    " max ", max(tr), "\n", sep = "")
cat("null  n=", length(nu), " min ", min(nu), " mean ", mean(nu),
    " max ", max(nu), "\n", sep = "")

set.seed(11)
B <- 20000L
for (k in c(2L, 3L, 4L, 5L)) {
  fail_mean <- 0L
  fail_min <- 0L
  for (b in seq_len(B)) {
    a <- sample(tr, k, replace = TRUE)
    z <- sample(nu, k, replace = TRUE)
    if (!(mean(a) > max(z))) fail_mean <- fail_mean + 1L
    if (!(min(a) > max(z))) fail_min <- fail_min + 1L
  }
  cat(sprintf(paste0("k=%d  mean(treat) > max(null) fails %d/%d;",
                     "  min(treat) > max(null) fails %d/%d\n"),
              k, fail_mean, B, fail_min, B))
}
## The margin the shipped seeds actually have is printed by the test
## itself when it fails, so this is only the rate at other seeds.
