# Reviewer, lane splinecurve, claim 1: the band at an unseen level.
#   Rscript splinecurve-rev-01-band.R > splinecurve-rev-log/01-band.txt
# Data: the curve-inference vignette's simulation, set.seed(4), 7200 rows
# over 20 subjects (same construction as dev/splinecurve-band.R).
LIB <- "C:/Users/adf44/source/r/wt-splinecurve-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
cat("frmtmb", format(packageVersion("frmtmb")), "frmtmb.spline",
    format(packageVersion("frmtmb.spline")), "from",
    dirname(find.package("frmtmb.spline")), "\n\n")

set.seed(4)
n_sub <- 20; n_rep <- 12; n_t <- 30
peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
sub <- rep(seq_len(n_sub), each = n_rep * n_t)
d <- data.frame(subject = factor(sub),
                t = rep(seq(0, 1, length.out = n_t), times = n_sub * n_rep))
h_sub <- rnorm(n_sub, 1, 0.12)
s_sub <- rnorm(n_sub, 0.5, 0.04)
d$v <- peak(d$t, h_sub[sub], s_sub[sub]) + rnorm(nrow(d), 0, 0.06)
d$f <- factor(ifelse(as.integer(d$subject) <= 10, "a", "b"))

lev_u <- c(levels(d$subject), "population")
grid <- data.frame(t = seq(0, 1, length.out = 80),
                   subject = factor("population", levels = lev_u))
seen <- transform(grid, subject = factor("3", levels = levels(d$subject)))
strip <- function(x) {
  attr(x, "fit") <- NULL
  x
}
try_msg <- function(expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) {
    paste0("REFUSED [", paste(class(r), collapse = "/"), "]: ",
           substr(conditionMessage(r), 1, 200))
  } else "ANSWERED"
}

## A. fs only: reproduce, and NA against NULL --------------------------------
fit <- frm(bf(v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)),
           family = gaussian(), data = d)
lb <- frm_lp_basis(fit, newdata = grid, re_formula = NA,
                   allow_new_levels = TRUE)
is_fs <- grepl("subject", lb$coef_names, fixed = TRUE)
cat("A. fs: ncol(A)", ncol(lb$A), " fs cols", sum(is_fs),
    " fs cols all exactly 0:", all(as.matrix(lb$A)[, is_fs] == 0),
    " extra_var all exactly 0:", all(lb$extra_var == 0), "\n")
cvNA <- frm_curve(fit, newdata = grid, re_formula = NA,
                  allow_new_levels = TRUE, nsim = 20000, seed = 1)
cvNULL <- frm_curve(fit, newdata = grid, re_formula = NULL,
                    allow_new_levels = TRUE, nsim = 20000, seed = 1)
cat("A. crit", format(cvNA$.crit_sim[1], digits = 10),
    " mcse", format(attr(cvNA, "check")$crit_mcse, digits = 6), "\n")
sNA <- strip(cvNA); sNULL <- strip(cvNULL)
attr(sNULL, "spec")$re_formula <- NA
cat("A. NA vs NULL at the unseen fs level: identical() after aligning the",
    "stored re_formula:", identical(sNA, sNULL), "\n")
A <- as.matrix(lb$A)
Ak <- A[, !is_fs, drop = FALSE]
se_k <- sqrt(diag(Ak %*% lb$V[!is_fs, !is_fs] %*% t(Ak)))
cat("A. .se identical() to hand restriction:",
    identical(cvNA$.se, unname(se_k)), "\n")
g <- mgcv::gam(v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5),
               data = d, method = "ML")
pg <- mgcv::predict.gam(g, newdata = transform(grid,
                          subject = factor("1", levels = levels(d$subject))),
                        exclude = "s(t,subject)", se.fit = TRUE)
r <- cvNA$.se / as.numeric(pg$se.fit)
cat("A. mgcv ML exclude=: max |est diff|",
    format(max(abs(cvNA$.estimate - pg$fit)), digits = 4),
    " se ratio range", format(range(r), digits = 5), " median",
    format(median(r), digits = 5), "\n")
# the population curve against the seen-subject curve's band width
cvS <- frm_curve(fit, newdata = seen, re_formula = NA, simultaneous = FALSE)
cat("A. median se: unseen", format(median(cvNA$.se), digits = 5),
    " seen subject 3", format(median(cvS$.se), digits = 5), "\n\n")

## B. fs AND (1 | subject) at an unseen subject ------------------------------
fit2 <- frm(bf(v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5) +
                 (1 | subject)),
            family = gaussian(), data = d)
vc <- VarCorr(fit2)
sd_int <- vc$subject$sd[1, 1]
cat("B. sd(subject intercept)", format(sd_int, digits = 8), "\n")
for (rf in list(NA, NULL)) {
  lb2 <- frm_lp_basis(fit2, newdata = grid, re_formula = rf,
                      allow_new_levels = TRUE)
  cv2 <- frm_curve(fit2, newdata = grid, re_formula = rf,
                   allow_new_levels = TRUE, simultaneous = FALSE)
  C2 <- as.matrix(lb2$A)
  coef_var <- diag(C2 %*% lb2$V %*% t(C2))
  cat("B. re_formula =", deparse(rf), ": extra_var range",
      format(range(lb2$extra_var), digits = 8),
      " max |extra_var / sd^2 - 1|",
      format(max(abs(lb2$extra_var / sd_int^2 - 1)), digits = 3),
      " max |.se^2 - coef_var - extra_var| / .se^2",
      format(max(abs(cv2$.se^2 - coef_var - lb2$extra_var) / cv2$.se^2),
             digits = 3),
      " median .se", format(median(cv2$.se), digits = 5), "\n")
}
cat("\n")

## C. s(subject, bs = "re") at an unseen level ---------------------------------
fit3 <- frm(bf(v ~ s(t, k = 12) + s(subject, bs = "re")),
            family = gaussian(), data = d)
cat("C. re smooth, unseen, TRUE, NA  :",
    try_msg(frm_curve(fit3, newdata = grid, re_formula = NA,
                      allow_new_levels = TRUE, simultaneous = FALSE)), "\n")
cat("C. re smooth, unseen, TRUE, NULL:",
    try_msg(frm_curve(fit3, newdata = grid, re_formula = NULL,
                      allow_new_levels = TRUE, simultaneous = FALSE)), "\n")
cat("C. re smooth, unseen, FALSE     :",
    try_msg(frm_curve(fit3, newdata = grid, simultaneous = FALSE)), "\n")
cat("C. re smooth, unseen, deriv TRUE:",
    try_msg(frm_curve_deriv(fit3, var = "t", newdata = grid,
                            allow_new_levels = TRUE, simultaneous = FALSE)),
    "\n")
cat("C. re smooth, SEEN, TRUE        :",
    try_msg(frm_curve(fit3, newdata = seen, allow_new_levels = TRUE,
                      simultaneous = FALSE)), "\n\n")

## D. t2(t, subject, bs = c("cr", "re")) at an unseen level --------------------
fit4 <- frm(bf(v ~ t2(t, subject, bs = c("cr", "re"), k = c(5, 20))),
            family = gaussian(), data = d)
cat("D. t2 re margin, unseen, TRUE   :",
    try_msg(frm_curve(fit4, newdata = grid, allow_new_levels = TRUE,
                      simultaneous = FALSE)), "\n")
cat("D. t2 re margin, unseen, FALSE  :",
    try_msg(frm_curve(fit4, newdata = grid, simultaneous = FALSE)), "\n")
cat("D. t2 re margin, feature TRUE   :",
    try_msg(frm_curve_feature(fit4, var = "t", type = "maximum",
                              newdata = grid, allow_new_levels = TRUE)), "\n")
cat("D. t2 re margin, SEEN, TRUE     :",
    try_msg(frm_curve(fit4, newdata = seen, allow_new_levels = TRUE,
                      simultaneous = FALSE)), "\n\n")

## E. fs with by = f -----------------------------------------------------------
fit5 <- frm(bf(v ~ f + s(t, subject, bs = "fs", by = f, k = 5)),
            family = gaussian(), data = d)
seen_f <- transform(seen, f = factor("a", levels = levels(d$f)))
grid_f <- transform(grid, f = factor("a", levels = levels(d$f)))
for (nm in c("seen, FALSE", "seen, TRUE", "unseen, TRUE", "unseen, FALSE")) {
  nd <- if (startsWith(nm, "seen")) seen_f else grid_f
  al <- endsWith(nm, "TRUE")
  for (rf in list(NA, NULL)) {
    cat("E. fs by=f,", nm, ", re_formula =", deparse(rf), ":",
        try_msg(frm_curve(fit5, newdata = nd, re_formula = rf,
                          allow_new_levels = al, simultaneous = FALSE)), "\n")
  }
}
cat("E. fs by=f, frm_linpred unseen TRUE:",
    try_msg(frm_linpred(fit5, newdata = grid_f, allow_new_levels = TRUE,
                        re_formula = NA)), "\n")
