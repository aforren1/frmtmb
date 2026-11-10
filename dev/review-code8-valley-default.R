# Review of lane code8: the DEFAULT frm() fit of y = A exp(B x), x on
# [x0, x0 + 1], against the best of several polishes of its own
# objective. What does each build say about the default stop?
#   RV_ARM=before|after Rscript dev/review-code8-valley-default.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
rows <- list()
for (x0 in c(30, 60, 100, 200)) for (seed in 1:8) {
  set.seed(seed)
  n <- 60
  d <- data.frame(x = runif(n, x0, x0 + 1))
  d$y <- 2 * exp((3 / x0) * d$x) + rnorm(n, 0, 1)
  r <- rfit(bf(y ~ A * exp(B * x), A ~ 1, B ~ 1, nl = TRUE), data = d)
  if (r$err) next
  f <- r$f
  obj <- f$obj
  x <- f$opt$par
  f0 <- obj$fn(x)
  best <- f0
  z <- x
  for (it in 1:4) {
    o1 <- nlminb(z, obj$fn, obj$gr, control = list(iter.max = 5000,
                                                    eval.max = 10000))
    if (o1$objective < best) { best <- o1$objective; z <- o1$par }
    o2 <- optim(z, obj$fn, obj$gr, method = "BFGS",
                control = list(maxit = 20000, reltol = 1e-15))
    if (o2$value < best) { best <- o2$value; z <- o2$par }
  }
  txt <- c(r$w, r$m)
  rows[[length(rows) + 1L]] <- data.frame(x0 = x0, seed = seed,
    code = f$opt$convergence, short = signif(f0 - best, 4),
    says_max = any(grepl("so the estimates are at the maximum", txt)),
    grad_w = any(grepl("Large maximum absolute gradient", txt)),
    at_max_m = any(r$mc == "frmtmb_converged_at_max"),
    nwarn = length(r$w), lost = length(se_of(f)$lost),
    conds = substr(cond_txt(r), 1, 100))
}
X <- do.call(rbind, rows)
options(width = 220)
print(X, right = FALSE)
write.table(X, sprintf("dev/review-code8-out/valley-default-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("short > 1e-3:", sum(X$short > 1e-3), " of those with no warning:",
    sum(X$short > 1e-3 & X$nwarn == 0), " saying at max:",
    sum(X$short > 1e-3 & (X$says_max | X$at_max_m)), "\n")
cat("DONE\n")
