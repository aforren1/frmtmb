# Lane splinecurve: what band does frm_curve() report for the POPULATION
# curve of v ~ s(t) + s(t, subject, bs = "fs") when the grid names an
# unseen subject with allow_new_levels = TRUE?
#
# Data: the curve-inference vignette's own simulation, set.seed(4).
# Run: Rscript dev/splinecurve-band.R > dev/splinecurve-band.log
lib <- "C:/Users/adf44/source/r/wt-splinecurve-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
cat("frmtmb.spline", format(packageVersion("frmtmb.spline")), "from",
    dirname(find.package("frmtmb.spline")), "\n\n")

set.seed(4)
n_sub <- 20
n_rep <- 12
n_t <- 30
peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
sub <- rep(seq_len(n_sub), each = n_rep * n_t)
d <- data.frame(
  subject = factor(sub),
  trial = rep(seq_len(n_sub * n_rep), each = n_t),
  t = rep(seq(0, 1, length.out = n_t), times = n_sub * n_rep))
h_sub <- rnorm(n_sub, 1, 0.12)
s_sub <- rnorm(n_sub, 0.5, 0.04)
d$v <- peak(d$t, h_sub[sub], s_sub[sub]) + rnorm(nrow(d), 0, 0.06)

fit <- frm(bf(v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)),
           family = gaussian(), data = d)

lev_u <- c(levels(d$subject), "population")
grid <- data.frame(t = seq(0, 1, length.out = 80),
                   subject = factor("population", levels = lev_u))

## 1. What core's seam returns at the unseen fs level --------------------
for (rf in list(NA, NULL)) {
  lb <- frm_lp_basis(fit, newdata = grid, re_formula = rf,
                     allow_new_levels = TRUE)
  cat("== frm_lp_basis(re_formula = ", deparse(rf),
      ", allow_new_levels = TRUE) at the unseen level\n", sep = "")
  cat("  columns of A:", ncol(lb$A), "\n")
  grp <- sub("[.][^.]*$", "", lb$coef_names)
  print(table(grp))
  cat("  extra_var: range", format(range(lb$extra_var)),
      " all exactly zero:", all(lb$extra_var == 0), "\n")
  is_fs <- grepl("subject", lb$coef_names, fixed = TRUE)
  cat("  fs columns of A:", sum(is_fs), " all exactly zero:",
      all(as.matrix(lb$A)[, is_fs] == 0), "\n\n")
}
lb <- frm_lp_basis(fit, newdata = grid, re_formula = NA,
                   allow_new_levels = TRUE)
cat("coef_names (first and last few):\n")
print(head(lb$coef_names, 16))
print(tail(lb$coef_names, 4))

## 2. frm_curve at the unseen level --------------------------------------
cv <- frm_curve(fit, newdata = grid, re_formula = NA, allow_new_levels = TRUE,
                nsim = 20000, seed = 1)
print(cv)
eta_ref <- as.numeric(frm_linpred(fit, newdata = grid, re_formula = NA,
                                  allow_new_levels = TRUE))
se_ref <- as.numeric(frm_linpred(fit, newdata = grid, re_formula = NA,
                                 allow_new_levels = TRUE,
                                 se.fit = TRUE)$se.fit)
cat("\nestimate identical() to frm_linpred(allow_new_levels = TRUE):",
    identical(cv$.estimate, eta_ref), " max abs diff",
    format(max(abs(cv$.estimate - eta_ref))), "\n")
cat("se vs frm_linpred(se.fit = TRUE): max |ratio - 1|",
    format(max(abs(cv$.se / se_ref - 1))), "\n")
cat("self-check cov_rel_error:", format(attr(cv, "check")$cov_rel_error),
    "\n")

## 3. (b) exact: the same seam restricted by hand to s(t) + intercept -----
A <- as.matrix(lb$A)
keep <- !grepl("subject", lb$coef_names, fixed = TRUE)
cat("\nkept columns:", sum(keep), "of", ncol(A), ":",
    paste(unique(sub("[.][^.]*$", "", lb$coef_names[keep])), collapse = ", "),
    "\n")
Ak <- A[, keep, drop = FALSE]
se_k <- sqrt(diag(Ak %*% lb$V[keep, keep] %*% t(Ak)))
cat("se restricted by hand vs frm_curve .se: max |ratio - 1|",
    format(max(abs(cv$.se / se_k - 1))), " identical():",
    identical(cv$.se, unname(se_k)), "\n")

## 4. (a) loose: mgcv's own population curve ------------------------------
if (requireNamespace("mgcv", quietly = TRUE)) {
  g <- mgcv::gam(v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5),
                 data = d, method = "ML")
  pg <- mgcv::predict.gam(g, newdata = transform(grid,
                            subject = factor("1", levels = levels(d$subject))),
                          exclude = "s(t,subject)", se.fit = TRUE)
  r_se <- cv$.se / as.numeric(pg$se.fit)
  cat("\nmgcv ML fit, predict.gam(exclude = \"s(t,subject)\", se.fit = TRUE)\n")
  cat("  estimate: max abs diff", format(max(abs(cv$.estimate - pg$fit))),
      " against the estimate's own range",
      format(diff(range(cv$.estimate))), "\n")
  cat("  se ratio frm/mgcv: range", format(range(r_se), digits = 4),
      " median", format(median(r_se), digits = 4), "\n")
  cat("  mgcv sp:", format(g$sp, digits = 4), "\n")
}

## 5. For contrast: a BAR term's unseen level does carry extra_var -------
fit_re <- frm(bf(v ~ s(t, k = 12) + (1 | subject)), family = gaussian(),
              data = d)
lb_re <- frm_lp_basis(fit_re, newdata = grid, re_formula = NULL,
                      allow_new_levels = TRUE)
vc <- VarCorr(fit_re)
cat("\n(1 | subject) at re_formula = NULL, unseen level: extra_var range",
    format(range(lb_re$extra_var)), "\n")
print(vc)

## 6. The difference guard: two DIFFERENT unseen levels of (1 | subject) --
gA <- transform(grid[1:5, ], subject = factor("A", levels = c(lev_u, "A", "B")))
gB <- transform(grid[1:5, ], subject = factor("B", levels = c(lev_u, "A", "B")))
a <- frmtmb.spline:::sp_one_basis(fit_re, gA, NULL, NULL, NULL, TRUE)
b <- frmtmb.spline:::sp_one_basis(fit_re, gB, NULL, NULL, NULL, TRUE)
cat("\nsp_same_latent() on two DIFFERENT unseen levels A and B:",
    frmtmb.spline:::sp_same_latent(fit_re, a, b),
    "\n  so without the guard the difference se would omit",
    "2 * extra_var =", format(2 * a$lb$extra_var[1]), "per row, against",
    "a coefficient part of", format(sum(((a$C - b$C) %*% a$lb$V) *
                                          (a$C - b$C))), "\n")
res <- tryCatch(frm_curve(fit_re, newdata = gA, contrast = gB,
                          re_formula = NULL, allow_new_levels = TRUE,
                          simultaneous = FALSE),
                error = function(e) conditionMessage(e))
cat("frm_curve(contrast = ) on it now:", if (is.character(res)) "REFUSED:" else
  "returned", if (is.character(res)) substr(res, 1, 120), "\n")
