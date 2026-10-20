# Re-check of lane illse2, finding 1: anova() between two whitened REML
# fits and profile intervals of a whitened REML fit, against lme4.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-rev-rcanova.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
set.seed(1)
deg <- 5
d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, length.out = 200)),
                h = factor(sample(1:10, 200, TRUE)))
d$y <- sin(3 * d$x) + rnorm(20, 0, 0.5)[d$g] + rnorm(10, 0, 0.3)[d$h] +
  rnorm(200, 0, 0.3)
for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
rhs <- paste0("x", 1:deg, collapse = " + ")
fo1 <- as.formula(paste("y ~", rhs, "+ (1 | g)"))
fo2 <- as.formula(paste("y ~", rhs, "+ (1 | g) + (1 | h)"))
f1 <- suppressWarnings(frm(fo1, data = d, REML = TRUE))
f2 <- suppressWarnings(frm(fo2, data = d, REML = TRUE))
l1 <- lme4::lmer(fo1, data = d, REML = TRUE)
l2 <- lme4::lmer(fo2, data = d, REML = TRUE)
cat(sprintf("whitened: %s %s | logLik frm %.8f %.8f | lme4 %.8f %.8f\n",
            !is.null(f1$beta_whiten), !is.null(f2$beta_whiten),
            as.numeric(logLik(f1)), as.numeric(logLik(f2)),
            as.numeric(logLik(l1)), as.numeric(logLik(l2))))
a <- cap(anova(f1, f2))
print(a$v); for (w in a$w) cat("W:", substr(w, 1, 200), "\n")
al <- anova(l1, l2, refit = FALSE)
print(al)
ci <- cap(confint(f2, method = "profile"))
if (is_err(ci$v)) cat("profile confint ERROR:", ci$v, "\n") else print(ci$v)
for (w in ci$w) cat("W:", substr(w, 1, 200), "\n")
cl <- cap(confint(l2, method = "profile", oldNames = FALSE))
print(cl$v)
