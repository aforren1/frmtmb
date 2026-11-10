# Lane code8 (copy of the review script dev/review-code8-ridge.R), item 4: the 41-draw ridge tally on a build, and a
# mutation that reintroduces the original defect (a straight step off a
# curved ridge taken as curvature: se_curvature_real() without its
# ratio test, so a one-sided loss above 2 grad_tol keeps the direction).
#   RV_ARM=after Rscript dev/review-code8-ridge.R <none|mutant> <label>
a <- commandArgs(TRUE)
variant <- a[1]; lab <- a[2]
source("dev/review-code8-common.R")
if (variant == "mutant") {
  m <- function(fit, p, free, dir_free, lambda) {
    tol <- fit$control$grad_tol %||% 1e-3
    obj <- fit$obj
    saved <- obj_state_save(obj)
    on.exit(obj_state_restore(obj, saved), add = TRUE)
    d <- numeric(length(p))
    d[free] <- dir_free * sqrt(4 * tol / lambda)
    f <- function(x) tryCatch(se_quiet(obj$fn(x)), error = function(e) NA_real_)
    f0 <- f(p)
    full <- f(p + d) + f(p - d) - 2 * f0
    is.finite(full) && full >= 2 * tol
  }
  environment(m) <- ns
  assignInNamespace("se_curvature_real", m, "frmtmb")
}
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
y0 <- dn$yn
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
rows <- list()
for (j in -20:20) {
  dn$yn <- y0 * (1 + 2e-15 * j)
  r <- rfit(form, data = dn)
  f <- r$f
  s <- se_of(f)
  rows[[length(rows) + 1L]] <- data.frame(
    j = j, obj = sprintf("%.8f", f$opt$objective), code = f$opt$convergence,
    a = signif(f$opt$par[3], 4), k = signif(f$opt$par[4], 4),
    se_a = signif(s$se[["a_(Intercept)"]], 4),
    se_c0 = signif(s$se[["c0_(Intercept)"]], 4),
    se_x = signif(s$se[["c0_x"]], 4),
    se_sigma = signif(s$se[["sigma_(Intercept)"]], 5),
    lost = paste(sort(sprintf("%s:%s", names(s$lost), s$lost)), collapse = " "),
    conds = substr(cond_txt(r), 1, 120), time = round(r$time, 2))
}
X <- do.call(rbind, rows)
write.table(X, sprintf("dev/code8-out/ridge-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
print(X[, c("j", "obj", "code", "a", "k", "se_a", "se_c0", "se_x",
            "se_sigma", "time")], right = FALSE)
cat("LOST SETS\n")
tb <- table(X$lost)
for (i in seq_along(tb)) cat(tb[[i]], ":", names(tb)[i], "\n")
cat("CONDS\n")
tb <- table(X$conds)
for (i in seq_along(tb)) cat(tb[[i]], ":", names(tb)[i], "\n")
cat("se_c0 range", range(X$se_c0), "se_sigma range", range(X$se_sigma), "\n")
cat("DONE\n")
