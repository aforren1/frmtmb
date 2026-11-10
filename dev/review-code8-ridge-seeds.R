# Review of lane code8, item 4: the ridge c0 + exp(a)^k beyond seed 77.
# The invariance (a0, k, theta_1) -> (s a0, k / s, theta_1 + log s) holds
# for any data, so every fit should lose exactly a, k and theta_1.
#   RV_ARM=before|after Rscript dev/review-code8-ridge-seeds.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
rows <- list()
for (seed in c(1:12, 77)) {
  set.seed(seed)
  dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
  dn$x <- rnorm(nrow(dn))
  u <- rnorm(8, 0, 0.6)[dn$g]
  dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
  y0 <- dn$yn
  for (j in c(0, -7, 11, 3)) {
    dn$yn <- y0 * (1 + 2e-15 * j)
    r <- rfit(form, data = dn)
    if (r$err) { cat(seed, j, "error\n"); next }
    s <- se_of(r$f)
    rows[[length(rows) + 1L]] <- data.frame(seed = seed, j = j,
      obj = sprintf("%.7f", r$f$opt$objective), code = r$f$opt$convergence,
      lost = paste(sort(names(s$lost)), collapse = ","),
      se_c0 = signif(s$se[["c0_(Intercept)"]], 3),
      se_sigma = signif(s$se[["sigma_(Intercept)"]], 4),
      grad_w = any(grepl("Large maximum", r$w)),
      conds = substr(cond_txt(r), 1, 80))
  }
}
X <- do.call(rbind, rows)
options(width = 220)
print(X, right = FALSE)
write.table(X, sprintf("dev/review-code8-out/ridge-seeds-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
want <- "a_(Intercept),k_(Intercept),theta_1"
cat("fits", nrow(X), "lost exactly a,k,theta_1:", sum(X$lost == want),
    " gradient warning:", sum(X$grad_w), "\n")
tb <- table(X$lost); for (i in seq_along(tb)) cat(tb[[i]], ":", names(tb)[i], "\n")
cat("DONE\n")
