# Lane gpselimit: the gpfield review's 30 gaussian approximate gp() fits
# (dev/reviews/gpfield-review-basis.R, same data), with each fit's
# conditions, the standard errors of sdgp and lscale and the lost reasons.
# Saves each row's estimates, SEs and conditions to
# dev/gpselimit-out/sweep-<arm><tag>-<row>.rds for a bitwise comparison.
# Usage: Rscript gpselimit-sweep.R <lane|base> [rows] [jitter k] [tag]
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]
rows <- if (length(a) >= 2 && nzchar(a[2])) eval(parse(text = a[2])) else 1:30
jit <- if (length(a) >= 3) as.numeric(a[3]) else 0
tag <- if (length(a) >= 4) a[4] else ""
gs_libs(arm)
for (row in rows) {
  s <- gs_sweep_row(row, jit)
  t0 <- proc.time()[["elapsed"]]
  q <- gs_fit(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian())
  fit <- q$fit
  if (inherits(fit, "error")) {
    cat("SWEEP", row, "ERROR", conditionMessage(fit), "\n"); next
  }
  tf <- proc.time()[["elapsed"]] - t0
  bk <- Filter(function(b) b$covstruct == "hsgp", fit$frame$re_blocks)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  sp <- bk$gp_spacing %||% NA
  t1 <- proc.time()[["elapsed"]]
  ss <- gs_ses(fit)
  ts <- proc.time()[["elapsed"]] - t1
  nm <- names(ss$se)
  tn <- which(nm == "theta" | grepl("^theta", nm))
  se_th <- ss$se[tn]
  saveRDS(list(row = row, est = fit$estimates, par = fit$opt$par,
               ll = as.numeric(logLik(fit)), se = ss$se, lost = ss$lost,
               wclass = q$wclass, warn = q$warn, mclass = q$mclass,
               msg = q$msg, sdr_conds = ss$conds),
          file.path(GS_OUT, sprintf("sweep-%s%s-%02d.rds", arm, tag, row)))
  cat(sprintf(paste0("SWEEP %2d %s n %3d ls0 %.1f k %2d | ll %.6f conv %d | ",
                     "ls/spacing %.3g sd %.4g | basis %d | se(theta) %s | ",
                     "lost [%s] | warn [%s] | msg [%s] | sdrcond %d | ",
                     "fit %.1fs se %.1fs\n"),
              row, arm, s$n, s$ls0, s$k, as.numeric(logLik(fit)),
              fit$opt$convergence, exp(th[2]) / sp, exp(th[1]),
              "frmtmb_gp_basis" %in% q$wclass,
              paste(format(se_th, digits = 4), collapse = ","),
              paste(paste0(names(ss$lost), "=", ss$lost), collapse = ","),
              paste(unique(q$wclass), collapse = ";"),
              paste(unique(q$mclass), collapse = ";"),
              length(ss$conds), tf, ts))
}
