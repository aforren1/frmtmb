# Punch round 1: summarise the recovery arms from the RDS files alone.
#
# Directories, in order of precedence for a seed that appears in more
# than one: recov5 (sequential reruns of failed seeds on the final
# build), recov4 (final build, phase3b-lib7), recov3 (phase3b-lib5).
# A seed's record is the first one WITHOUT an error in that order; every
# error seen on the way is counted and classified. Each record carries
# the path of the frmtmb.eam it was fitted with, and the build of each
# arm is printed from it.
#
# Usage: Rscript dev/phase3b-summarise2.R <arm> [<arm> ...]
# Output: dev/phase3b-log/summary2-<arm>.txt
arms <- commandArgs(trailingOnly = TRUE)
dirs <- c("dev/phase3b-log/recov5", "dev/phase3b-log/recov4",
          "dev/phase3b-log/recov3")

wilson <- function(k, n) {
  z <- 1.959964; p <- k / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n),
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n))
}
row_of <- function(name, truth, est, lwr, upr) {
  ok <- is.finite(est) & is.finite(lwr) & is.finite(upr) & is.finite(truth)
  k <- sum(lwr[ok] < truth[ok] & truth[ok] < upr[ok]); n <- sum(ok)
  w <- wilson(k, n)
  sprintf("%-22s %9.4f %9.4f %8.4f %4d/%-4d %6.1f%%  [%4.1f, %4.1f]",
          name, mean(truth[ok]), mean(est[ok]), sd(est[ok]) / sqrt(n), k, n,
          100 * k / n, 100 * w[1], 100 * w[2])
}
header <- sprintf("%-22s %9s %9s %8s %9s %7s  %s", "quantity", "truth",
                  "mean", "mcse", "covered", "rate", "Wilson 95")

for (arm in arms) {
  files <- unlist(lapply(dirs, list.files,
                         pattern = sprintf("^%s-[0-9]+[.]rds$", arm),
                         full.names = TRUE))
  seedof <- function(f) as.integer(sub(".*-([0-9]+)[.]rds$", "\\1", f))
  seeds <- sort(unique(seedof(files)))
  seeds <- seeds[seeds < 9000]
  recs <- list(); errs <- character(0); builds <- character(0)
  for (s in seeds) {
    got <- NULL
    for (dd in dirs) {
      f <- file.path(dd, sprintf("%s-%d.rds", arm, s))
      if (!file.exists(f)) next
      r <- readRDS(f)
      if (!is.null(r$error)) {
        errs <- c(errs, if (grepl("bad_alloc", r$error)) "bad_alloc" else
          if (grepl("NA/NaN gradient", r$error)) "NaN gradient" else
            substr(r$error, 1, 40))
        next
      }
      got <- r; builds <- c(builds, basename(dirname(r$eam_path %||% "?/?")))
      break
    }
    if (!is.null(got)) recs[[length(recs) + 1L]] <- got
  }
  out <- c(sprintf("arm %s: %d seeds (%d to %d), %d with a fit", arm,
                   length(seeds), min(seeds), max(seeds), length(recs)),
           sprintf("errors met on the way (every attempt): %s",
                   if (length(errs)) paste(names(table(errs)), table(errs),
                                           sep = " ", collapse = "; ") else "none"),
           sprintf("builds of the fits used: %s",
                   paste(names(table(builds)), table(builds), sep = " ",
                         collapse = "; ")))
  x <- recs
  conv <- vapply(x, function(r) r$conv, 0)
  flag <- vapply(x, function(r) {
    length(r$diagnose) > 3L &&
      !any(grepl("No convergence problems detected", r$diagnose))
  }, TRUE)
  out <- c(out, sprintf("convergence code 0: %d of %d; diagnose() reported a finding: %d",
                        sum(conv == 0), length(x), sum(flag)))
  why <- unlist(lapply(x, function(r) {
    d <- r$diagnose[-(1:3)]
    unique(substr(d[nzchar(trimws(d)) & !grepl("^ ", d)], 1, 48))
  }))
  if (length(why)) {
    tb <- sort(table(why), decreasing = TRUE)
    out <- c(out, sprintf("  %3d  %s", as.integer(tb), names(tb)))
  }
  get <- function(row, col) vapply(x, function(r) {
    ci <- r$confint
    if (is.matrix(ci) && row %in% rownames(ci)) ci[row, col] else NA_real_
  }, 0)
  ub <- vapply(x, function(r) r$ndt_ub, 0)
  ndt_truth <- suppressWarnings(qlogis(0.25 / ub))
  ndt_truth[0.25 >= ub] <- NA
  n <- length(x)
  out <- c(out, "", header,
           row_of("mu intercept", rep(0.4, n), get("(Intercept)", "est"),
                  get("(Intercept)", "lwr"), get("(Intercept)", "upr")),
           row_of("mu condition", rep(0.9, n), get("cond", "est"),
                  get("cond", "lwr"), get("cond", "upr")),
           row_of("log bs", rep(log(1.4), n), get("bs_(Intercept)", "est"),
                  get("bs_(Intercept)", "lwr"), get("bs_(Intercept)", "upr")),
           row_of("ndt (link scale)", ndt_truth, get("ndt_(Intercept)", "est"),
                  get("ndt_(Intercept)", "lwr"), get("ndt_(Intercept)", "upr")),
           sprintf("ndt natural scale: mean %.4f (truth 0.25); truth below the bound on %d of %d",
                   mean(plogis(get("ndt_(Intercept)", "est")) * ub, na.rm = TRUE),
                   sum(is.finite(ndt_truth)), n),
           row_of("log sd(mu | s)", rep(log(0.35), n), get("theta_1", "est"),
                  get("theta_1", "lwr"), get("theta_1", "upr")),
           row_of("log sd(log bs | s)", rep(log(0.20), n), get("theta_2", "est"),
                  get("theta_2", "lwr"), get("theta_2", "upr")))
  if (arm %in% c("cont", "contfix", "contdl")) {
    le <- get("lambda_(Intercept)", "est")
    nc <- vapply(x, function(r) r$n_cont / r$n, 0)
    out <- c(out,
             row_of("lambda (logit), vs 0.05", rep(qlogis(0.05), n), le,
                    get("lambda_(Intercept)", "lwr"),
                    get("lambda_(Intercept)", "upr")),
             row_of("lambda (logit), vs share", qlogis(nc), le,
                    get("lambda_(Intercept)", "lwr"),
                    get("lambda_(Intercept)", "upr")),
             sprintf("lambda natural: mean %.4f; realized share of recorded rows: mean %.4f",
                     mean(plogis(le), na.rm = TRUE), mean(nc)),
             sprintf("window used: top mean %.3f, bottom mean %.4f",
                     mean(vapply(x, function(r) r$crange[2], 0)),
                     mean(vapply(x, function(r) r$crange[1], 0))))
  }
  if (arm == "collapse") {
    le <- get("lambda_(Intercept)", "est")
    dll <- vapply(x, function(r) r$loglik - (r$loglik_plain %||% NA), 0)
    edge <- le < -15
    nm <- vapply(x, function(r) any(grepl("end of its link: lambda",
                                           r$diagnose)), TRUE)
    out <- c(out, sprintf("lambda at the edge (logit below -15): %d of %d; diagnose() names it: %d; of those at the edge: %d",
                          sum(edge, na.rm = TRUE), n, sum(nm), sum(nm & edge, na.rm = TRUE)),
             sprintf("interior lambda: %s",
                     if (any(!edge, na.rm = TRUE)) paste(format(range(plogis(le[!edge])), digits = 3), collapse = " to ") else "none"),
             sprintf("|logLik - plain logLik| at the edge: max %.3g over %d with a plain fit",
                     max(abs(dll[edge]), na.rm = TRUE), sum(edge & is.finite(dll), na.rm = TRUE)))
  }
  if (arm == "cleft") {
    out <- c(out, sprintf("left-censored share mean %.4f; interval share mean %.4f",
                          mean(vapply(x, function(r) r$n_left / r$n, 0)),
                          mean(vapply(x, function(r) r$n_int / r$n, 0))))
  }
  writeLines(out, sprintf("dev/phase3b-log/summary2-%s.txt", arm))
  cat(out, sep = "\n"); cat("\n")
}
