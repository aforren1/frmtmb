# Item 5222: hmm(2, skew_normal()) with both states started at alpha = 0
# (dev/silent-hmm-probe.R, ported). The escape is one pass, so one alpha
# stays at 0. Also mixture() of two skew normals started the same way.
#   Rscript dev/quadseam2-hmm.R <libpath|base> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.latent)})
cat("lib:", dirname(find.package("frmtmb")), "\n")
cat("latent:", dirname(find.package("frmtmb.latent")), "\n")
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:3
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
ea <- function(f) {
  b <- f$estimates$betad
  signif(b[grep("^alpha", names(b))], 4)
}
fo <- bf(y ~ 1)
for (seed in seeds) {
  dd <- mk(seed)
  hf <- hmm(2, skew_normal(), time = t, group = ID)
  u <- frm(fo, family = hf, data = dd, dry_run = "objective")
  free <- suppressWarnings(frm(fo, family = hf, data = dd))
  bd <- u$estimates$betad
  bd[grep("^alpha", names(bd))] <- 0
  t0 <- proc.time()[["elapsed"]]
  stuck <- suppressWarnings(frm(fo, family = hf, data = dd,
                                start = list(betad = bd)))
  t1 <- proc.time()[["elapsed"]] - t0
  cat(sprintf(paste0("HMM seed %d free %.4f alpha %s | stall %.4f alpha %s ",
                     "escape %s secs %.1f\n"),
              seed, as.numeric(logLik(free)), paste(ea(free), collapse = " "),
              as.numeric(logLik(stuck)), paste(ea(stuck), collapse = " "),
              paste(format(stuck$opt$stationary_escape, digits = 4),
                    collapse = " "), t1))
}
# mixture(skew_normal, skew_normal), both alphas started at 0
for (seed in seeds) {
  set.seed(100 + seed)
  n <- 400
  z <- rbinom(n, 1, 0.5)
  y <- ifelse(z == 1, sim_sn(n, -3, 1, 5), sim_sn(n, 3, 1, -5))
  dd <- data.frame(y = y)
  mf <- mixture(skew_normal(), skew_normal())
  fo2 <- bf(y ~ 1)
  u <- frm(fo2, family = mf, data = dd, dry_run = "objective")
  free <- suppressWarnings(frm(fo2, family = mf, data = dd))
  bd <- u$estimates$betad
  bd[grep("^alpha", names(bd))] <- 0
  stuck <- suppressWarnings(frm(fo2, family = mf, data = dd,
                                start = list(betad = bd)))
  cat(sprintf(paste0("MIX seed %d free %.4f alpha %s | stall %.4f alpha %s ",
                     "escape %s\n"),
              seed, as.numeric(logLik(free)), paste(ea(free), collapse = " "),
              as.numeric(logLik(stuck)), paste(ea(stuck), collapse = " "),
              paste(format(stuck$opt$stationary_escape, digits = 4),
                    collapse = " ")))
}
cat("END\n")
