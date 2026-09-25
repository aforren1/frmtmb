# Two follow-ups from rv-features: why a hurdle draw at the fitted rows
# is not identical() though numerically equal, and whether se() at
# newdata sets the residual spread from newdata's se column.
source("dev/simnewdata-review/rv-prelude.R")
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n))
eta <- 0.3 + 0.5 * d$x
d$yh <- ifelse(runif(n) < 0.3, 0, 1 + rpois(n, exp(eta)))
fh <- frm(bf(yh ~ x, hu ~ x), family = hurdle_poisson(), data = d)
a <- simulate(fh, nsim = 3, seed = 4)
b <- simulate(fh, nsim = 3, seed = 4, newdata = d)
cat("hurdle identical:", identical(a, b), " all.equal:",
    isTRUE(all.equal(a, b)), "\n")
cat("  types:", typeof(a[[1]]), typeof(b[[1]]), "\n")
cat("  attributes a:", paste(names(attributes(a)), collapse = ","),
    " b:", paste(names(attributes(b)), collapse = ","), "\n")
cat("  rownames equal:", identical(rownames(a), rownames(b)), "\n")
print(all.equal(a, b))
d$se <- runif(n, 0.2, 1)
d$ym <- eta + rnorm(n, 0, d$se)
for (sg in c(FALSE, TRUE)) {
  fm <- frm(bf(ym | se(se, sigma = sg) ~ x), data = d)
  nd <- data.frame(x = 0, se = c(0.1, 2))
  s <- as.matrix(simulate(fm, nsim = 20000, seed = 1, newdata = nd))
  sig <- if (sg) sigma(fm) else 0
  cat(sprintf("se(sigma = %s): draw sd %.3f %.3f; expected %.3f %.3f\n",
              sg, sd(s[1, ]), sd(s[2, ]), sqrt(0.1^2 + sig^2),
              sqrt(4 + sig^2)))
}
