# Lane selimit2: the limit rule's internals on zi ~ x seed 5.
#   Rscript dev/selimit2-dbg.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(5)
n <- 200
d <- data.frame(x = stats::rnorm(n), N = sample(5:15, n, TRUE))
d$y <- stats::rbinom(n, d$N, stats::plogis(0.3 + 0.5 * d$x))
f <- suppressWarnings(frm(bf(y | trials(N) ~ x, zi ~ x),
                          family = zero_inflated_binomial(), data = d))
p <- f$opt$par
H <- f$obj$he(p)
om <- ns$outer_par_map(f)
cand <- ns$se_limit_candidates(f, om)
u <- rep_len(f$par_units %||% 1, length(p))
u[cand] <- ns$se_limit_units(f, om)[cand]
cat("u:", u, "\n")
cat("diag H nat:", abs(diag(H)) * u^2, "\n")
jd <- ns$se_limit_joint_dirs(f, H, p, u, cand, om)
cat("joint dirs:", length(jd), "\n")
for (w in jd) {
  cat("  w:", signif(w, 4), " z:", attr(w, "z"), "\n")
  f0 <- f$obj$fn(p)
  for (t in c(-12, -2, 2, 14)) {
    cat(sprintf("    t %+g: %.4g\n", t, f$obj$fn(p + t * w) - f0))
  }
}
for (j in 3:4) {
  f0 <- f$obj$fn(p)
  for (t in c(-12, -2, 2, 14)) {
    q <- p
    q[j] <- p[j] + t * u[j]
    cat(sprintf("  coord %d t %+g: %.4g\n", j, t, f$obj$fn(q) - f0))
  }
}
el <- ns$se_edge_limit(f, H)
print(el)
str(ns$se_limit_info(f))
