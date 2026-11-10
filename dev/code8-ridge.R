# Part C: the SE verdict on test-se-tier.R's ridge (c0 + exp(a)^k, a ~ 1 +
# (1 | g), seed 77), refit with y * (1 + 2e-15 * j), j = -20..20.
#   CODE8_ARM=before|after Rscript dev/code8-ridge.R <label> [trace]
a <- commandArgs(TRUE)
lab <- a[1]
do_trace <- length(a) > 1 && a[2] == "trace"
source("dev/code8-common.R")
options(width = 200)
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
y0 <- dn$yn
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
if (do_trace) {
  ns <- asNamespace("frmtmb")
  trace("se_ridge_curvature", exit = quote(cat("  ridge_curv k=", k,
    " ev=", format(e$values[k], digits = 3), " ->",
    format(returnValue(), digits = 3), "\n")), where = ns, print = FALSE)
  trace("se_curvature_real", exit = quote(cat("  curv_real lambda=",
    format(lambda, digits = 3), " ->", returnValue(), "\n")), where = ns,
    print = FALSE)
}
rows <- list()
for (j in -20:20) {
  dn$yn <- y0 * (1 + 2e-15 * j)
  if (do_trace) cat("== draw", j, "\n")
  r <- rfit(form, data = dn)
  f <- r$f
  nm <- fns$outer_par_names(f)
  p <- f$opt$par
  sdr <- quiet(fns$sdr_of(f))
  se <- sqrt(diag(as.matrix(sdr$cov.fixed)))
  lost <- sdr$se_lost
  rows[[length(rows) + 1L]] <- data.frame(
    j = j, obj = f$opt$objective, code = f$opt$convergence,
    a = p[match("a_(Intercept)", nm)], k = p[match("k_(Intercept)", nm)],
    theta = p[match("theta_1", nm)],
    se_a = se[match("a_(Intercept)", nm)],
    se_sigma = se[match("sigma", nm)],
    lost = paste(sprintf("%s:%s", names(lost), lost), collapse = " "),
    conds = substr(cond_txt(r), 1, 160))
}
X <- do.call(rbind, rows)
write.table(X, sprintf("dev/code8-out/ridge-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
X$obj <- sprintf("%.8f", X$obj)
print(X[, c("j", "obj", "code", "a", "k", "theta", "se_a", "se_sigma",
            "lost")], right = FALSE, digits = 4)
cat("lost sets:\n"); print(table(X$lost))
cat("se_a finite:", sum(is.finite(X$se_a)), "of", nrow(X), "\n")
cat("conds:\n"); print(table(X$conds))
cat("DONE\n")
