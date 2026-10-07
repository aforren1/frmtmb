# Reviewer of lane setier: the degree-7 raw polynomial with (1 | g)
# (seed 7, x on [1, 2]): which SEs are kept, against lmer() and against
# the exact-Hessian fit without (1 | g); the tier 3 spectrum and tau.
#   Rscript dev/setier-rev-deg7.R lane|base [deg]
args <- commandArgs(TRUE)
arm <- args[1]
deg <- if (length(args) > 1) as.integer(args[2]) else 7L
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
set.seed(deg)
d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
rhs <- paste0("x", 1:deg, collapse = " + ")
f <- suppressWarnings(frm(as.formula(paste("y ~", rhs, "+ (1 | g)")), data = d))
fe <- suppressWarnings(fixef(f))
m <- lmer(as.formula(paste("y ~", rhs, "+ (1 | g)")), data = d, REML = FALSE)
print(summary(m)$optinfo$warnings)
cat("lmer fixef cols:", names(fixef(m)), "\n")
sl <- sqrt(diag(as.matrix(vcov(m))))
cat(sprintf("%-12s %12s %12s %12s %12s\n", "coef", "frmtmb", "SE", "lmer est",
            "lmer SE"))
for (i in seq_len(nrow(fe))) {
  cat(sprintf("%-12s %12.5g %12.5g %12.5g %12.5g\n", rownames(fe)[i],
              fe[i, 1], fe[i, 2], fixef(m)[i], sl[i]))
}
# the same data with the group means as known offsets: exact Hessian ref
lm0 <- lm(as.formula(paste("y ~", rhs)), data = d)
cat("lm SE (no g):", signif(sqrt(diag(vcov(lm0))), 4), "\n")
H <- f$cache$hessian_fixed$H
if (is.null(H)) H <- ns$fit_outer_hessian(f)$H
E <- f$cache$hessian_fixed$E
an <- ns$cov_from_hessian(f, H, E)
cat("tier", an$tier, "lost", names(an$lost), "\n")
cat("ev (unit diag):", signif(an$ev, 3), "\n")
cat("cond-SE of kept vs solve(H) SE:\n")
print(signif(rbind(kept = sqrt(diag(an$V)),
                   solveH = suppressWarnings(sqrt(diag(solve(H))))), 4))
