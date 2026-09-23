# Lane wt-reunc: inspect the random-effect block structure and the
# base behavior of fitted() Est.Error at a known level.
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/rellib-r3"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(1)
G <- 8; m <- 5
d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                h = factor(rep(seq_len(4), times = G * m / 4)),
                x = rnorm(G * m))
d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 0.7)[d$g] + rnorm(4, 0, 0.5)[d$h] +
  rnorm(G * m, 0, 1)
fit <- frm(bf(y ~ x + (1 + x | g) + (1 | h)), data = d)
for (bk in fit$frame$re_blocks) {
  cat("block covstruct", bk$covstruct, "dim", bk$dim, "nlev",
      length(bk$levels), "\n")
  cat("  b_idx", head(bk$b_idx), " c_idx", head(bk$c_idx), "\n")
  for (cp in bk$components) {
    cat("  comp label", cp$label, "lp_key", cp$lp_key, "offset", cp$offset,
        "dim", cp$dim, "cnms", cp$cnms, "bar", deparse(cp$bar), "\n")
  }
}
lp <- fit$frame$linpreds[[1]]
cat("Z dim", dim(lp$Z), class(lp$Z), "\n")
sdr <- frmtmb:::sdr_of(fit)
Q <- sdr$jointPrecision
cat("Q null?", is.null(Q), "\n")
if (!is.null(Q)) print(table(rownames(Q)))
print(frmtmb:::outer_par_map(fit))
print(names(fit$estimates))
print(fit$REML)
# fitted() Est.Error at a known level versus the base prediction
f1 <- fitted(fit)
cat("fitted Est.Error first 3:", f1[1:3, 2], "\n")
p0 <- frm_linpred(fit, se.fit = TRUE)
cat("se.fit first 3:", p0$se.fit[1:3], "\n")
pn <- frm_linpred(fit, se.fit = TRUE, re_formula = NA)
cat("se.fit NA first 3:", pn$se.fit[1:3], "\n")
