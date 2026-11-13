# Does hmm(K, skew_normal()) build, and does a fit started on alpha = 0
# stay there? Arm "base" reads rellib-r20; arm "lane" reads the lane lib.
# Usage: Rscript dev/silent-hmm-probe.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
}
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.latent)})
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")),
    "from", find.package("frmtmb"), "\n")
cat("latent from", find.package("frmtmb.latent"), "\n")

fam <- hmm(2, skew_normal())
cat("post$stationary names:",
    paste(names(fam$post$stationary), collapse = ","), "\n")

sim_sn <- function(n, mu, sigma, alpha) {
  d <- alpha / sqrt(1 + alpha^2)
  z0 <- abs(rnorm(n)); z1 <- rnorm(n)
  mu + sigma * (d * z0 + sqrt(1 - d^2) * z1)
}
mk <- function(seed, nseq = 20, T = 40, alpha = 5) {
  set.seed(seed)
  P <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE)
  out <- vector("list", nseq)
  for (s in seq_len(nseq)) {
    st <- integer(T); st[1] <- sample(1:2, 1)
    for (t in 2:T) st[t] <- sample(1:2, 1, prob = P[st[t - 1], ])
    y <- ifelse(st == 1, sim_sn(T, -3, 1, alpha), sim_sn(T, 3, 1, -alpha))
    out[[s]] <- data.frame(ID = s, t = seq_len(T), y = y)
  }
  do.call(rbind, out)
}
fo <- bf(y ~ 1)
for (seed in 1:3) {
  dd <- mk(seed)
  hf <- hmm(2, skew_normal(), time = t, group = ID)
  u <- frm(fo, family = hf, data = dd, dry_run = "objective")
  if (seed == 1) {
    cat("start names:\n"); print(u$estimates)
  }
  free <- frm(fo, family = hf, data = dd)
  bd <- u$estimates$betad
  bd[grep("^alpha", names(bd))] <- 0
  stuck <- frm(fo, family = hf, data = dd, start = list(betad = bd))
  ea <- function(f) {
    b <- f$estimates$betad
    signif(b[grep("^alpha", names(b))], 4)
  }
  cat(sprintf("seed %d free logLik %.4f alpha %s | from alpha=0 logLik %.4f alpha %s escape %s\n",
              seed, as.numeric(logLik(free)),
              paste(ea(free), collapse = " "),
              as.numeric(logLik(stuck)),
              paste(ea(stuck), collapse = " "),
              paste(format(stuck$opt$stationary_escape), collapse = " ")))
}
