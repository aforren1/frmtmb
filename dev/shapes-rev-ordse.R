# Reviewer: the ordinal interop standard error, diagnosed.
#
# marginaleffects reports 0.00029 for the MIDDLE category's average
# marginal effect where MASS::polr reports 0.02016 on the same data with
# the same point estimate. This asks whether the cause is the one the
# shape work exposes: the thresholds are not in the coefficient vector
# the interop seams hand out, so the numeric jacobian has no threshold
# direction and the delta method loses the variance that lives there.
#
# The check: redo the delta method HERE over the whole outer parameter
# vector, using vcov(fit, full = TRUE), which does carry the
# thresholds. If that lands on polr's number, the diagnosis holds.

.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
source(file.path(TREE, "dev/shapes-rev-fixtures.R"))
dd <- rev_data()
f <- frm(bf(ord ~ x) + cumulative(), data = dd)

cat("estimated_coef_names(): ",
    paste(frmtmb:::estimated_coef_names(f), collapse = ", "), "\n")
cat("vcov_estimated() dim:   ",
    paste(dim(frmtmb::vcov_estimated(f)), collapse = " x "), "\n")
cat("vcov(full = TRUE) rows: ",
    paste(rownames(suppressWarnings(vcov(f, full = TRUE))),
          collapse = ", "), "\n")
cat("insight::get_parameters(): ",
    paste(insight::get_parameters(f)$Parameter, collapse = ", "), "\n")
cat("marginaleffects::get_coef(): ",
    paste(names(marginaleffects::get_coef(f)), collapse = ", "), "\n\n")

# the AME of x for each category, by a central difference over the
# whole outer parameter vector
ame <- function(fit) {
  p <- fitted(fit)                       # n x 4 x K
  est <- p[, "Estimate", , drop = TRUE]
  h <- 1e-5
  d2 <- fit$frame[["data"]]
  nd <- transform(model.frame(fit), x = model.frame(fit)$x + h)
  colMeans((fitted(fit, newdata = nd)[, "Estimate", ] - est) / h)
}
V <- suppressWarnings(vcov(f, full = TRUE))
v0 <- frmtmb:::fit_outer_vector(f)
map <- frmtmb:::outer_par_map(f)
a0 <- ame(f)
J <- matrix(NA_real_, length(a0), length(v0))
for (j in seq_along(v0)) {
  hh <- 1e-5 * max(1, abs(v0[j]))
  vp <- v0; vp[j] <- vp[j] + hh
  vm <- v0; vm[j] <- vm[j] - hh
  J[, j] <- (ame(frmtmb:::fit_set_outer(f, vp, map)) -
               ame(frmtmb:::fit_set_outer(f, vm, map))) / (2 * hh)
}
se_full <- sqrt(diag(J %*% V %*% t(J)))
# and the same with ONLY the coefficient block, which is what the
# interop seams hand marginaleffects
keep <- which(rownames(V) %in% frmtmb:::estimated_coef_names(f))
se_beta <- sqrt(diag(J[, keep, drop = FALSE] %*%
                       V[keep, keep, drop = FALSE] %*%
                       t(J[, keep, drop = FALSE])))
cat("category   AME        se(full outer vector)   se(coef block only)\n")
for (k in seq_along(a0))
  cat(sprintf("  %-8d %9.5f  %14.5f %20.5f\n", k, a0[k], se_full[k],
              se_beta[k]))
cat("\nMASS::polr for comparison\n")
pl <- MASS::polr(ord ~ x, data = dd, Hess = TRUE)
s <- as.data.frame(suppressWarnings(marginaleffects::avg_slopes(pl)))
print(s[, c("group", "estimate", "std.error")])
