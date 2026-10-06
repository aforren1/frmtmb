# Reviewer: find a converged fit with random effects whose lost standard
# error is a FIXED effect (separation), test-se-check.R's separation
# design plus a random intercept. Usage:
# Rscript dev/cifixrev-sepre.R <lib or "base"> <seed> [<seed> ...]
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
f6 <- function(v) format(signif(v, 6))
warns <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}
ctl <- frmtmb_control(optCtrl = list(eval.max = 4000, iter.max = 4000))
for (s in as.integer(args[-1])) {
  set.seed(s)
  n <- 240
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(rep(1:12, each = 20)),
                  f = factor(rep(c("a", "b"), 120)))
  d$yb <- as.integer(d$z > 0)
  # half the outcomes are not separated, so the intercept, x and the
  # random intercept stay identified
  flip <- d$f == "b"
  d$zb <- ifelse(flip, 0, d$z)
  d$yb[flip] <- stats::rbinom(sum(flip), 1,
                              stats::plogis(0.5 * d$x[flip] +
                                              stats::rnorm(12)[d$g[flip]]))
  r <- warns(frm(yb ~ x + f + zb + (1 | g), family = bernoulli(), data = d,
                 control = ctl))
  fit <- r$value
  sdr <- ns$sdr_of(fit)
  cat("\n== seed", s, "code", fit$opt$convergence, "evals", fit$opt$evals,
      "lost:", paste(names(sdr$se_lost), sdr$se_lost), "\n")
  for (w in r$warnings) cat("  W:", substr(w, 1, 200), "\n")
  if (!length(sdr$se_lost)) next
  jc <- ns$get_joint_cov(fit)
  Q <- ns$joint_precision(fit)
  cat("lost_pos", jc$lost_pos, rownames(Q)[jc$lost_pos], "\n")
  keep <- setdiff(seq_len(nrow(Q)), jc$lost_pos)
  if (length(jc$lost_pos)) {
    Vind <- as.matrix(Matrix::solve(Q[keep, keep]))
    cat("repaired vs solve(Q[-lost, -lost]):",
        f6(max(abs(jc$V[keep, keep] - Vind)) / max(abs(Vind))), "\n")
  }
  nd <- data.frame(x = c(0, 1, 0, 0), f = factor(c("a", "a", "b", "b"),
                                                 levels = c("a", "b")),
                   zb = c(0, 0, 0, 1), g = factor(c(1, 2, 3, 4)))
  for (rf in list(NULL, NA)) {
    pw <- warns(frm_linpred(fit, newdata = nd, se.fit = TRUE,
                            re_formula = rf))
    cat("frm_linpred re_formula =", deparse(rf), "se:",
        f6(pw$value$se.fit), "warnings:", length(pw$warnings), "\n")
    for (w in pw$warnings) cat("  W:", substr(w, 1, 160), "\n")
  }
  pp <- warns(predict(fit, newdata = nd))
  cat("predict Est.Error:", f6(pp$value[, "Est.Error"]), "warnings:",
      length(pp$warnings), "\n")
  for (w in pp$warnings) cat("  W:", substr(w, 1, 160), "\n")
  rc <- warns(ranef(fit, condVar = TRUE))
  cv <- as.data.frame(rc$value)
  bpos <- which(rownames(Q) == "b")
  cat("ranef condsd / sqrt(diag V_bb) range:",
      f6(range(cv$condsd / sqrt(diag(jc$V)[bpos]))), "warnings:",
      length(rc$warnings), "\n")
  vc <- warns(VarCorr(fit))
  print(as.data.frame(vc$value))
  cat("VarCorr warnings:", length(vc$warnings), "\n")
  ce <- warns(conditional_effects(fit, "x", re_formula = NULL))
  cat("conditional_effects(x, re_formula = NULL) se__ NaN:",
      sum(is.nan(ce$value[[1]]$se__)), "of", nrow(ce$value[[1]]),
      "warnings:", length(ce$warnings), "\n")
  for (w in ce$warnings) cat("  W:", substr(w, 1, 160), "\n")
  ce2 <- warns(conditional_effects(fit, "zb", re_formula = NA))
  cat("conditional_effects(zb, re_formula = NA) se__ NaN:",
      sum(is.nan(ce2$value[[1]]$se__)), "of", nrow(ce2$value[[1]]),
      "warnings:", length(ce2$warnings), "\n")
  em <- warns(as.data.frame(emmeans::emmeans(fit, ~ f)))
  cat("emmeans SE:", f6(em$value$SE), "warnings:", length(em$warnings),
      "\n")
  for (w in em$warnings) cat("  W:", substr(w, 1, 160), "\n")
  print(summary(fit)$fixed)
  jcp <- frm_joint_cov(fit)
  cat("frm_joint_cov NaN rows:", which(apply(is.nan(jcp$V), 1, all)), "\n")
}
