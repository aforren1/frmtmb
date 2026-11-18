# Review of lane sparsese: follow-ups on the "base" design of
# sparsese-review-designs.R (degree-7 raw polynomial + f + (1 | g), REML,
# 41 joint parameters), dense (limit 1e6) and sparse (limit 20) arms.
#  1. emmeans(~ f) SE against frm_linpred(se.fit) at emmeans's grid;
#  2. predict() draws: Est.Error against sqrt(se.fit^2 + sigma^2) (seen
#     levels) and with a new level, 20000 draws;
#  3. frm_curve() on each arm, and on the sparse arm the error it gives.
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(1)
ng <- 30L
n <- 300L
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)),
                h = factor(sample(rep(seq_len(12), length.out = n))),
                f = factor(rep(c("a", "b"), each = 1, length.out = n)))
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
ug <- stats::rnorm(ng, 0, 0.5)
invisible(stats::rnorm(ng, 0, 0.4))
invisible(stats::rnorm(12, 0, 0.3))
d$y <- sin(3 * d$x) + ug[d$g] + 0.3 * (d$f == "b") + stats::rnorm(n, 0, 0.3)
fo <- y ~ f + x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
lims <- c(dense = 1e6, sparse = 20)
for (arm in names(lims)) {
  options(frmtmb.joint_root_max = lims[[arm]])
  fit <- q(frm(fo, data = d, REML = TRUE))
  invisible(q(frmtmb:::sdr_of(fit)))
  cat("== arm", arm, " rows attr:",
      length(attr(frmtmb:::sdr_of(fit)$joint_root, "rows")), "\n")
  em <- q(as.data.frame(emmeans::emmeans(fit, ~ f)))
  grid <- data.frame(f = factor(c("a", "b"), levels = c("a", "b")))
  for (k in 1:7) grid[[paste0("x", k)]] <- mean(d[[paste0("x", k)]])
  lp <- q(frm_linpred(fit, newdata = grid, re_formula = NA, se.fit = TRUE))
  cat("emmeans SE:", format(em$SE, digits = 10), "\n")
  cat("linpred SE:", format(lp$se.fit, digits = 10), "\n")
  cat("emmeans / linpred - 1:", format(em$SE / lp$se.fit - 1, digits = 3),
      " estimate diff:", format(max(abs(em$emmean - lp$fit)), digits = 3),
      "\n")
  nd <- d[c(1:6, 151:156), ]
  sig <- exp(fit$opt$par[["betad"]])
  se <- q(frm_linpred(fit, newdata = nd, se.fit = TRUE))$se.fit
  set.seed(5)
  pr <- q(predict(fit, newdata = nd, ndraws = 20000))
  cat("predict Est.Error / sqrt(se^2 + sigma^2) - 1:",
      format(range(pr[, "Est.Error"] / sqrt(se^2 + sig^2) - 1), digits = 3),
      "\n")
  sdg <- exp(fit$opt$par[["theta"]])
  ndn <- nd
  ndn$g <- factor(rep("new", nrow(nd)))
  sen <- q(frm_linpred(fit, newdata = ndn, se.fit = TRUE,
                       allow_new_levels = TRUE))$se.fit
  se0 <- q(frm_linpred(fit, newdata = ndn, se.fit = TRUE, re_formula = NA))$se.fit
  cat("linpred new level se^2 - (noRE se^2 + sd^2), rel:",
      format(max(abs(sen^2 - se0^2 - sdg^2) / sen^2), digits = 3), "\n")
  set.seed(6)
  prn <- q(predict(fit, newdata = ndn, ndraws = 20000, allow_new_levels = TRUE,
                   sample_new_levels = "gaussian"))
  cat("predict new level Est.Error / sqrt(se0^2 + sd^2 + sigma^2) - 1:",
      format(range(prn[, "Est.Error"] / sqrt(se0^2 + sdg^2 + sig^2) - 1),
             digits = 3), "\n")
  cv <- tryCatch(q(frm_curve(fit, newdata = grid[1, , drop = FALSE][rep(1, 5), ],
                             seed = 1)),
                 error = function(e) e)
  if (inherits(cv, "error")) {
    cat("frm_curve ERROR:", substr(conditionMessage(cv), 1, 400), "\n")
  } else {
    cat("frm_curve ok; cov_rel_error",
        format(attr(cv, "check")$cov_rel_error, digits = 3), "\n")
  }
  lb <- q(frm_lp_basis(fit, newdata = nd, re_formula = NA))
  v1 <- rowSums((lb$A %*% lb$V) * lb$A)
  se_na <- q(frm_linpred(fit, newdata = nd, re_formula = NA, se.fit = TRUE))$se.fit
  cat("frm_lp_basis: V_root", !is.null(lb$V_root),
      "; sqrt(A V A') / se.fit - 1:",
      format(max(abs(sqrt(pmax(v1, 0)) / se_na - 1)), digits = 3), "\n")
  rm(fit)
}
