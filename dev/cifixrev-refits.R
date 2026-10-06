# Reviewer: joint_cov_repair() on fits with random effects whose lost
# directions are not coordinate axes or are fixed effects.
# Usage: Rscript dev/cifixrev-refits.R <lib or "base">
# R1: test-se-check.R's gr(g, by = f) fixture (seed 11): theta_2 and
#     theta_3 lost on a NEGATIVE eigen-direction that x loads 0.10 on.
# R2: the nonlinear ridge y ~ a + b, a ~ 0 + f, b ~ 1 with (1 | g) on b:
#     a lost FIXED-EFFECT direction shared by several coefficients.
# R3: ls ~ mo(income) * age + (1 | g), se_mo_data(71) and (12).
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
se_mo_data <- function(seed) {
  set.seed(seed)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d$g <- factor(rep(1:10, 10))
  d$ls <- d$ls + rnorm(10, 0, 4)[d$g]
  d
}
report <- function(tag, fw, nd, dpar = NULL, ce = NULL, em = NULL) {
  fit <- fw$value
  sdr <- ns$sdr_of(fit)
  cat("\n==", tag, "code", fit$opt$convergence, "fit warnings",
      length(fw$warnings), "\n")
  for (w in fw$warnings) cat("  W:", substr(w, 1, 200), "\n")
  cat(tag, "lost:", paste(names(sdr$se_lost), sdr$se_lost), "\n")
  nul <- sdr$se_null
  if (NCOL(nul)) {
    cat(tag, "null columns, nonzero entries per column:",
        colSums(nul != 0), "\n")
  }
  Q <- ns$joint_precision(fit)
  jc <- ns$get_joint_cov(fit)
  cat(tag, "lost_pos:", jc$lost_pos, rownames(Q)[jc$lost_pos], "\n")
  if (length(jc$lost_pos)) {
    keep <- setdiff(seq_len(nrow(Q)), jc$lost_pos)
    Vind <- as.matrix(Matrix::solve(Q[keep, keep]))
    cat(tag, "repaired vs solve(Q[-lost, -lost]): max abs diff / max:",
        f6(max(abs(jc$V[keep, keep] - Vind)) / max(abs(Vind))), "\n")
    rr <- diag(jc$V)[keep] / diag(Vind)
    cat(tag, "diag ratio repaired / conditional range:", f6(range(rr)),
        "\n")
    ev <- eigen(jc$V, symmetric = TRUE, only.values = TRUE)$values
    cat(tag, "repaired eigen min/max:", f6(min(ev) / max(ev)), "\n")
  }
  ev0 <- eigen(as.matrix(Matrix::solve(Q)), symmetric = TRUE,
               only.values = TRUE)$values
  cat(tag, "solve(Q) eigen min/max:", f6(min(ev0) / max(ev0)), "\n")
  for (rf in list(NULL, NA)) {
    pw <- warns(frm_linpred(fit, newdata = nd, se.fit = TRUE,
                            re_formula = rf, dpar = dpar))
    cat(tag, "frm_linpred re_formula =", deparse(rf), "se:",
        f6(pw$value$se.fit), "warnings:", length(pw$warnings), "\n")
    for (w in pw$warnings) cat("  W:", substr(w, 1, 120), "\n")
  }
  rc <- warns(ranef(fit, condVar = TRUE))
  cv <- as.data.frame(rc$value)
  bpos <- which(rownames(Q) == "b")
  if (length(bpos) == nrow(cv)) {
    cat(tag, "ranef condsd / sqrt(diag V_bb) range:",
        f6(range(cv$condsd / sqrt(diag(jc$V)[bpos]))), "\n")
  }
  vc <- warns(VarCorr(fit))
  cat(tag, "VarCorr warnings:", length(vc$warnings), "\n")
  if (!is.null(ce)) {
    cw <- warns(conditional_effects(fit, ce, re_formula = NULL))
    s <- cw$value[[1]]$se__
    cat(tag, "conditional_effects(", ce, ", re_formula = NULL) se__ range:",
        f6(range(s)), "NaN:", sum(is.nan(s)), "of", length(s),
        "warnings:", length(cw$warnings), "\n")
    for (w in cw$warnings) cat("  W:", substr(w, 1, 120), "\n")
  }
  if (!is.null(em)) {
    ew <- warns(as.data.frame(emmeans::emmeans(fit, em)))
    cat(tag, "emmeans SE:", f6(ew$value$SE), "warnings:",
        length(ew$warnings), "\n")
    for (w in ew$warnings) cat("  W:", substr(w, 1, 120), "\n")
  }
}

# R1
set.seed(11)
dd <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, 10)))
dd$f <- factor(ifelse(as.integer(dd$g) <= 8, "a", "b"))
u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
dd$y <- stats::rnorm(160, 1 + 0.5 * dd$x + u[dd$g, 1] +
                       u[dd$g, 2] * dd$x, 1)
f1 <- warns(frm(bf(y ~ x + (1 + x | gr(g, by = f))), family = gaussian(),
                data = dd))
report("R1", f1, dd[c(1, 2, 9, 10), ], ce = "x")

# R2
set.seed(1)
k <- 10
d2 <- data.frame(f = factor(rep(seq_len(k), each = 6)),
                 g = factor(rep(1:6, k)))
d2$y <- rnorm(k)[d2$f] + rnorm(6, 0, 0.5)[d2$g] + rnorm(k * 6, 0, 0.5)
f2 <- warns(frm(bf(y ~ a + b, a ~ 0 + f, b ~ 1 + (1 | g), nl = TRUE),
                data = d2))
report("R2", f2, d2[c(1, 7), ], dpar = "a")

# R3
for (s in c(71, 12, 7)) {
  f3 <- warns(frm(ls ~ mo(income) * age + (1 | g), data = se_mo_data(s)))
  report(paste0("R3s", s), f3, se_mo_data(s)[1:4, ], ce = "income",
         em = ~ income)
}
