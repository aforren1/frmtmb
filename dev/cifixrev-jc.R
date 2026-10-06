# Reviewer: attack joint_cov_repair() (lane cifix).
# Usage: Rscript dev/cifixrev-jc.R <lib or "base">
# F1: the test fixture (seed 21, n = 100), two smoothing sds lost.
#     Repaired V against an independent inverse of Q with the lost rows
#     removed; the Schur complement of Q against the outer Hessian.
# F2: a GLMM whose covariate separates the outcome: random effects AND
#     a lost fixed-effect direction. Predictions, warnings, ranef condVar,
#     VarCorr, conditional_effects(), emmeans.
# F3: mgcv on a smooth whose smoothing parameter goes to the boundary.
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
g <- function(nm) get(nm, envir = ns)
f6 <- function(v) format(signif(v, 6))
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
warns <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

## F1 -------------------------------------------------------------------
set.seed(21)
n <- 100
d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                fac = factor(rep(c("A", "B"), length.out = n)))
d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
  stats::rnorm(n, 0, 0.3)
fit <- suppressWarnings(
  frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)), data = d))
sdr <- g("sdr_of")(fit)
cat("\n== F1 lost:", paste(names(sdr$se_lost), sdr$se_lost), "\n")
cat("F1 theta:", f6(fit$estimates$theta), "\n")
Q <- g("joint_precision")(fit)
jc <- g("get_joint_cov")(fit)
p <- nrow(Q)
r <- fit$obj$env$random
fo <- setdiff(seq_len(p), r)
cat("F1 lost_pos:", jc$lost_pos, " names:", rownames(Q)[jc$lost_pos], "\n")
nul <- sdr$se_null
cat("F1 se_null columns are coordinate vectors:",
    all(colSums(nul != 0) == 1), " ncol:", NCOL(nul), "\n")
if (length(jc$lost_pos)) {
  keep <- setdiff(seq_len(p), jc$lost_pos)
  Vind <- as.matrix(Matrix::solve(Q[keep, keep]))
  cat("F1 repaired vs solve(Q[-lost, -lost]): max abs diff / max abs:",
      f6(rel(jc$V[keep, keep], Vind)), "\n")
  dr <- diag(jc$V)[keep] / diag(Vind) - 1
  cat("F1 diag ratio - 1 range:", f6(range(dr)), "\n")
  cat("F1 repaired V lost rows all zero:",
      all(jc$V[jc$lost_pos, ] == 0), "\n")
  ev <- eigen(jc$V, symmetric = TRUE, only.values = TRUE)$values
  cat("F1 repaired V eigen min / max:", f6(min(ev) / max(ev)), "\n")
}
ev0 <- eigen(as.matrix(Matrix::solve(Q)), symmetric = TRUE,
             only.values = TRUE)$values
cat("F1 solve(Q) eigen range:", f6(range(ev0)), "\n")
# Schur complement of Qrr: the H that sdreport() was handed
S <- as.matrix(Q[fo, fo] - Q[fo, r] %*% Matrix::solve(Q[r, r], Q[r, fo]))
H <- g("fit_outer_hessian")(fit)$H
cat("F1 Schur(Q) vs fit_outer_hessian H: max abs diff / max abs:",
    f6(rel(S, H)), "\n")
gx <- d$x[-1] - diff(d$x) / 2
nd <- data.frame(x = rep(gx, 2),
                 fac = factor(rep(c("A", "B"), each = length(gx))))
lp <- warns(frm_linpred(fit, newdata = nd, se.fit = TRUE))
cat("F1 frm_linpred se range:", f6(range(lp$value$se.fit)),
    " warnings:", length(lp$warnings), "\n")
pr <- warns(frm_linpred(fit, newdata = nd, se.fit = TRUE, re_formula = NA))
cat("F1 frm_linpred(re_formula = NA) se range:", f6(range(pr$value$se.fit)),
    " warnings:", length(pr$warnings), "\n")
ce <- warns(conditional_effects(fit, "x"))
cat("F1 conditional_effects se__ range:",
    f6(range(ce$value[[1]]$se__)), " warnings:", length(ce$warnings), "\n")
rcv <- warns(ranef(fit, condVar = TRUE))
cat("F1 ranef condVar warnings:", length(rcv$warnings), "\n")

## F2 -------------------------------------------------------------------
set.seed(5)
ng <- 20
n2 <- 400
d2 <- data.frame(g = factor(rep(seq_len(ng), each = n2 / ng)),
                 z = stats::rnorm(n2),
                 f = factor(sample(c("a", "b", "c"), n2, TRUE)))
d2$x <- stats::rnorm(n2)
u <- stats::rnorm(ng, 0, 0.8)[d2$g]
eta <- -0.3 + 0.7 * d2$z + c(a = 0, b = 0.4, c = -0.5)[d2$f] + u
d2$y <- stats::rbinom(n2, 1, stats::plogis(eta))
# x separates within f == "c": every "c" row with x > 0 is a success
d2$x[d2$f == "c"] <- abs(d2$x[d2$f == "c"]) *
  ifelse(d2$y[d2$f == "c"] == 1, 1, -1)
d2$xc <- ifelse(d2$f == "c", d2$x, 0)
fit2w <- warns(frm(bf(y ~ z + f + xc + (1 | g)), family = bernoulli(),
                   control = frmtmb_control(optCtrl = list(eval.max = 4000,
                                                  iter.max = 4000)),
                   data = d2))
fit2 <- fit2w$value
cat("\n== F2 fit warnings:", length(fit2w$warnings), "\n")
for (w in fit2w$warnings) cat("  W:", substr(w, 1, 160), "\n")
sdr2 <- g("sdr_of")(fit2)
cat("F2 convergence:", fit2$opt$convergence, " lost:",
    paste(names(sdr2$se_lost), sdr2$se_lost), "\n")
cat("F2 beta:", f6(fit2$estimates$beta), "\n")
jc2 <- g("get_joint_cov")(fit2)
cat("F2 lost_pos:", jc2$lost_pos, " null cols:", NCOL(jc2$null), "\n")
Q2 <- g("joint_precision")(fit2)
if (length(jc2$lost_pos)) {
  keep <- setdiff(seq_len(nrow(Q2)), jc2$lost_pos)
  Vind <- as.matrix(Matrix::solve(Q2[keep, keep]))
  cat("F2 repaired vs solve(Q[-lost, -lost]):",
      f6(rel(jc2$V[keep, keep], Vind)), "\n")
}
nd2 <- data.frame(z = c(0, 0, 0, 0), f = factor(c("a", "b", "c", "c"),
                                               levels = c("a", "b", "c")),
                  xc = c(0, 0, 0, 1), g = factor(c(1, 2, 3, 4)))
for (rf in list(NULL, NA)) {
  pw <- warns(frm_linpred(fit2, newdata = nd2, se.fit = TRUE,
                      re_formula = rf))
  cat("F2 frm_linpred re_formula =", deparse(rf), " se:",
      f6(pw$value$se.fit), " warnings:", length(pw$warnings), "\n")
  for (w in pw$warnings) cat("  W:", substr(w, 1, 160), "\n")
}
lw <- warns(frm_linpred(fit2, newdata = nd2, se.fit = TRUE))
cat("F2 frm_linpred se:", f6(lw$value$se.fit), " warnings:",
    length(lw$warnings), "\n")
rc <- warns(ranef(fit2, condVar = TRUE))
cv <- as.data.frame(rc$value)
cat("F2 ranef condsd range:", f6(range(cv$condsd)), " warnings:",
    length(rc$warnings), "\n")
bpos <- which(rownames(Q2) == "b")
cat("F2 sqrt(diag(V_bb)) of get_joint_cov range:",
    f6(range(sqrt(diag(jc2$V)[bpos]))), "\n")
cat("F2 condsd / sqrt(diag V_bb) range:",
    f6(range(cv$condsd / sqrt(diag(jc2$V)[bpos]))), "\n")
vc <- warns(VarCorr(fit2))
print(as.data.frame(vc$value)); cat("F2 VarCorr warnings:",
    length(vc$warnings), "\n")
cew <- warns(conditional_effects(fit2, "z", re_formula = NULL))
cat("F2 conditional_effects(re_formula = NULL) se__ range:",
    f6(range(cew$value[[1]]$se__)), " NaN:", sum(is.nan(cew$value[[1]]$se__)),
    " warnings:", length(cew$warnings), "\n")
for (w in cew$warnings) cat("  W:", substr(w, 1, 160), "\n")
if (requireNamespace("emmeans", quietly = TRUE)) {
  em <- warns(as.data.frame(emmeans::emmeans(fit2, ~ f)))
  cat("F2 emmeans SE:", f6(em$value$SE), " warnings:",
      length(em$warnings), "\n")
  for (w in em$warnings) cat("  W:", substr(w, 1, 160), "\n")
}

## F3 -------------------------------------------------------------------
set.seed(3)
n3 <- 200
d3 <- data.frame(x = stats::runif(n3), z = stats::runif(n3))
d3$y <- 1 + 2 * d3$x + sin(2 * pi * d3$z) + stats::rnorm(n3, 0, 0.3)
fit3w <- warns(frm(bf(y ~ s(x) + s(z)), data = d3))
fit3 <- fit3w$value
sdr3 <- g("sdr_of")(fit3)
cat("\n== F3 theta:", f6(fit3$estimates$theta), " lost:",
    paste(names(sdr3$se_lost), sdr3$se_lost), " fit warnings:",
    length(fit3w$warnings), "\n")
nd3 <- data.frame(x = seq(0.05, 0.95, length.out = 7), z = 0.5)
s3 <- warns(frm_linpred(fit3, newdata = nd3, se.fit = TRUE))$value$se.fit
cat("F3 frmtmb se:", f6(s3), "\n")
if (requireNamespace("mgcv", quietly = TRUE)) {
  gm <- mgcv::gam(y ~ s(x) + s(z), data = d3, method = "ML")
  cat("F3 mgcv sp:", f6(gm$sp), "\n")
  Xp <- mgcv::predict.gam(gm, nd3, type = "lpmatrix")
  sp_ <- sqrt(rowSums((Xp %*% gm$Vp) * Xp))
  sc_ <- sqrt(rowSums((Xp %*% gm$Vc) * Xp))
  cat("F3 mgcv Vp se:", f6(sp_), "\n")
  cat("F3 mgcv Vc se:", f6(sc_), "\n")
  cat("F3 frmtmb / mgcv Vp:", f6(s3 / sp_), "\n")
}
