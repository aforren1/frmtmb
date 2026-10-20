# Re-check of lane illse2, finding 1: a profile interval of a whitened
# REML fit, against lme4's.
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
set.seed(1)
deg <- 5
d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, length.out = 200)),
                h = factor(sample(1:10, 200, TRUE)))
d$y <- sin(3 * d$x) + rnorm(20, 0, 0.5)[d$g] + rnorm(10, 0, 0.3)[d$h] +
  rnorm(200, 0, 0.3)
for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
fo2 <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + "),
                        "+ (1 | g) + (1 | h)"))
f2 <- suppressWarnings(frm(fo2, data = d, REML = TRUE))
w <- suppressWarnings(confint(f2))
print(rownames(w))
for (pm in rownames(w)) {
  r <- cap(confint(f2, parm = pm, method = "profile"))
  if (is_err(r$v)) cat(pm, "ERROR:", r$v, "\n") else print(r$v)
  for (x in r$w) cat("  W:", substr(x, 1, 160), "\n")
}
l2 <- lme4::lmer(fo2, data = d, REML = TRUE)
cat("lme4 REML sigma", sigma(l2), "sds", attr(lme4::VarCorr(l2)$g, "stddev"),
    attr(lme4::VarCorr(l2)$h, "stddev"), "\n")
print(exp(confint(f2)))
