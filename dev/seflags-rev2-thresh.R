# Reviewer of lane seflags, re-check of B2: the "half its curvature"
# rule of se_span_remeasure() on identified raw polynomials (degrees 6
# to 8 on [1, 2], with and without (1 | g), seeds 1 to 20). Per fit:
# code, lost reasons, max |SE / ref - 1| over the kept coefficients
# against lm() / lmer(), and every re-measured direction's remaining
# share (c0 - mix) / c0 (the rule keeps it above 0.5). The function is
# the lane's own with one recording line added.
#   Rscript dev/seflags-rev2-thresh.R <lane|op2> <out.tsv>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev2-common.R")
out <- commandArgs(trailingOnly = TRUE)[2]
rev_span_log <- list()
txt <- deparse(ns$se_span_remeasure)
i <- grep("c0 - mix > se_span_keep * c0", txt, fixed = TRUE)
stopifnot(length(i) == 1L)
txt <- append(txt, "rev_span_log[[length(rev_span_log) + 1L]] <<- c(c0, mix)",
              after = i - 1L)
f2 <- eval(parse(text = txt))
environment(f2) <- ns
utils::assignInNamespace("se_span_remeasure", f2, ns = "frmtmb")
cat("deg\tre\tseed\tcode\tlost\tn_lost\tmax_rel_kept\tshares\n", file = out)
for (deg in 6:8) for (re in c(FALSE, TRUE)) for (s in 1:20) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:deg, collapse = " + ")
  fo <- stats::as.formula(paste("y ~", rhs, if (re) "+ (1 | g)"))
  cn <- c("(Intercept)", paste0("x", 1:deg))
  ref <- if (re) {
    m <- suppressWarnings(suppressMessages(lme4::lmer(fo, data = d,
                                                      REML = FALSE)))
    sqrt(diag(as.matrix(stats::vcov(m))))
  } else sqrt(diag(stats::vcov(stats::lm(fo, data = d))))
  rev_span_log <- list()
  f <- suppressWarnings(frm(fo, data = d))
  l <- suppressWarnings(ns$sdr_of(f))$se_lost
  se <- se_outer(f)[cn]
  kept <- is.finite(se)
  rel <- if (any(kept)) max(abs(se[kept] / ref[kept] - 1)) else NA
  sh <- vapply(rev_span_log, function(v) (v[1] - v[2]) / v[1], 0)
  cat(deg, re, s, f$opt$convergence,
      if (length(l)) paste(names(table(l)), table(l), collapse = ",") else "-",
      length(l), signif(rel, 4),
      if (length(sh)) paste(signif(sh, 3), collapse = ",") else "-",
      sep = "\t", file = out, append = TRUE)
  cat("\n", file = out, append = TRUE)
}
cat("done\n")
