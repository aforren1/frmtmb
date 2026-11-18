# rel084-kt copy of dev/reviews/gpselimit-short-a.R: sources dev/rel084-kt/gpselimit-review-lib.R (rellib-r23); OUT = dev/rel084-log/kt/gpselimit-out. No other change.
# Short check of lane gpselimit's fix round 1, items 1 to 3, on one arm.
#  S  the 30 sweep rows, default start (13 white, 17 resolved)
#  A  the plateau stops: rows 6, 12, 18, 24, 30 started at 3e-5 of the
#     spacing
#  Y  gp(x, by = f) on 5 white rows (seeds 9600 + row), j = 0
#  N  the review's nl2 case, b0 -0.8, k 20, seeds 2 and 4
# Per fit: conditions (type, first class, gp text), lost reasons, SEs;
# saved to gpselimit-short-out/a-<arm>.rds for gpselimit-short-cmp.R.
# Usage: Rscript gpselimit-short-a.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/gpselimit-review-lib.R")
arm <- commandArgs(TRUE)[1]
rv_libs(arm)
ns <- asNamespace("frmtmb")
OUT <- file.path(RV_WT, "dev/rel084-log/kt/gpselimit-out")
dir.create(OUT, showWarnings = FALSE)
res <- list()
gp_conds <- function(conds) {
  g <- Filter(function(c) "frmtmb_gp_basis" %in% c$class, conds)
  if (!length(g)) return("none")
  paste(vapply(g, function(c) {
    paste0(substr(c$type, 1, 1), ":",
           if (grepl("stopped short", c$text)) "short" else
             if (grepl("^Boundary fit", c$text)) "boundary-msg" else "plain")
  }, ""), collapse = ";")
}
keep <- function(lab, q) {
  fit <- q$fit
  if (inherits(fit, "error")) {
    cat("FIT", lab, arm, "ERROR", conditionMessage(fit), "\n")
    return(invisible(NULL))
  }
  ss <- rv_ses(fit)
  wh <- vapply(rv_hsgp(fit), function(bk) {
    ns$gp_basis_white(bk, fit$estimates$theta[bk$theta_idx])
  }, NA)
  cat(sprintf(paste0("FIT %-16s %s | ll %.6f conv %d | white %s | ",
                     "gp [%s] | lost [%s] | conds [%s]\n"),
              lab, arm, as.numeric(logLik(fit)), fit$opt$convergence,
              paste(as.integer(wh), collapse = ""), gp_conds(q$conds),
              paste(paste0(names(ss$lost), "=", ss$lost), collapse = ","),
              rv_cstr(q$conds)))
  res[[lab]] <<- list(est = fit$estimates, par = fit$opt$par,
                      ll = as.numeric(logLik(fit)), se = ss$se,
                      lost = ss$lost, conds = q$conds, white = wh)
  invisible(fit)
}
for (row in 1:30) {
  s <- rv_sweep_row(row)
  keep(sprintf("S%02d", row),
       rv_fit(frm(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian())))
}
for (row in c(6, 12, 18, 24, 30)) {
  s <- rv_sweep_row(row)
  f0 <- rv_q(frm(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian()))
  sp <- rv_hsgp(f0)[[1]]$gp_spacing
  keep(sprintf("A%02d-3e-5", row),
       rv_fit(frm(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian(),
                  start = list(theta = c(log(30), log(3e-5 * sp))))))
}
for (row in c(1, 3, 9, 21, 25)) {
  s <- rv_sweep_row(row)
  s$d$f <- factor(rep(c("a", "b"), length.out = nrow(s$d)))
  fit <- keep(sprintf("Y%02d", row),
              rv_fit(frm(y ~ gp(x, by = f, k = s$k), data = s$d)))
  if (is.null(fit)) next
  r <- res[[sprintf("Y%02d", row)]]
  thp <- which(names(r$se) == "theta")
  wh <- unlist(lapply(rv_hsgp(fit), function(bk) {
    if (ns$gp_basis_white(bk, fit$estimates$theta[bk$theta_idx])) {
      thp[bk$theta_idx]
    }
  }))
  oth <- setdiff(seq_along(r$se), wh)
  cat(sprintf(paste0("   BY %02d %s white theta %d, NaN %d | other SEs %d, ",
                     "finite %d\n"), row, arm, length(wh),
              sum(is.nan(r$se[wh])), length(oth),
              sum(is.finite(r$se[oth]))))
}
for (seed in c(2, 4)) {
  set.seed(52300 + seed)
  n <- 120L
  x <- sort(round(stats::runif(n, 0, 10), 3))
  z <- round(stats::rnorm(n), 3)
  f <- 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  z <- abs(z) * 2
  y <- (2 + fs) * (1 - exp(0.8 * z)) + 0.1 * stats::rnorm(n)
  d <- data.frame(x, z, y)
  keep(sprintf("N-seed%d", seed),
       rv_fit(frm(bf(y ~ a * (1 - exp(-b * z)), a ~ gp(x, k = 20), b ~ 1,
                     nl = TRUE), data = d, family = gaussian())))
  cl <- vapply(res[[sprintf("N-seed%d", seed)]]$conds, function(c) {
    paste0(substr(c$type, 1, 1), ":", c$class[1])
  }, "")
  cat("   order", paste(cl, collapse = " > "), "\n")
}
saveRDS(res, file.path(OUT, sprintf("a-%s.rds", arm)))
