# Review of lane illse2, item 3: sd_up_escape() where the true sd is 0.
# Ordinary designs, gaussian and poisson with (1 | g) and no group
# effect, against lme4 (ML / Laplace); lane against 0.71.0.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-rev-esc.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
out_rds <- sprintf("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-rev-log/esc-%s.rds", arm)
rows <- list()
one <- function(fam, seed) {
  set.seed(seed)
  n <- 200
  d <- data.frame(x = rnorm(n), g = factor(rep(1:20, length.out = n)))
  if (fam == "gaussian") {
    d$y <- 1 + 0.5 * d$x + rnorm(n)
    t0 <- proc.time()[[3]]
    r <- cap(frm(y ~ x + (1 | g), data = d))
    tt <- proc.time()[[3]] - t0
    l <- cap(lme4::lmer(y ~ x + (1 | g), data = d, REML = FALSE))
  } else {
    d$y <- rpois(n, exp(0.3 + 0.4 * d$x))
    t0 <- proc.time()[[3]]
    r <- cap(frm(bf(y ~ x + (1 | g)) + poisson(), data = d))
    tt <- proc.time()[[3]] - t0
    l <- cap(lme4::glmer(y ~ x + (1 | g), data = d, family = poisson))
  }
  f <- r$v
  p <- f$opt$par
  se <- suppressWarnings(sqrt(diag(vcov(f))))
  data.frame(fam = fam, seed = seed, code = f$opt$convergence,
             ll = as.numeric(logLik(f)), ll_l = as.numeric(logLik(l$v)),
             sd = exp(p[names(p) == "theta"])[1],
             sd_l = attr(lme4::VarCorr(l$v)$g, "stddev")[[1]],
             esc = if (is.null(f$opt$sd_escape)) NA else f$opt$sd_escape[[1]],
             w = wtag(r$w), lost = lost_of(f),
             se = paste(format(se, digits = 10), collapse = " "),
             time = tt)
}
for (fam in c("gaussian", "poisson")) {
  for (s in 1:40) rows[[length(rows) + 1L]] <- one(fam, s)
}
tab <- do.call(rbind, rows)
print(tab[, c("fam", "seed", "code", "sd", "sd_l", "esc", "w", "lost")])
for (fm in unique(tab$fam)) {
  s <- tab[tab$fam == fm, ]
  cat(sprintf("%s (%s): n %d | sd < 1e-3: %d | escapes %d | warnings %d | lost %d | ll - lme4: min %.3e max %.3e | time %.2f s total\n",
              fm, arm, nrow(s), sum(s$sd < 1e-3), sum(!is.na(s$esc)),
              sum(s$w != "-"), sum(s$lost != "none"),
              min(s$ll - s$ll_l), max(s$ll - s$ll_l), sum(s$time)))
}
saveRDS(tab, out_rds)
