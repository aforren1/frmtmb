# REVIEW of lane eamhier, item 2.1. An INDEPENDENT reader of the
# per-replicate records in dev/eamhier-rec/, written so that no number
# below is taken from the lane's summary scripts.
#
# Run: Rscript --vanilla dev/eamhierrev-records.R

rev_read <- function(dir) {
  fs <- sort(list.files(dir, pattern = "[.]tsv$", full.names = TRUE))
  rows <- lapply(fs, function(f) {
    ln <- readLines(f, warn = FALSE)
    ln <- ln[nzchar(ln)]
    kv <- strsplit(ln[[1L]], "\t", fixed = TRUE)[[1L]]
    out <- as.list(trimws(sub("^[^=]*=", "", kv)))
    names(out) <- sub("=.*$", "", kv)
    out$file <- basename(f)
    out$mtime <- format(file.info(f)$mtime, "%Y-%m-%d %H:%M:%S")
    out
  })
  nms <- unique(unlist(lapply(rows, names)))
  d <- do.call(rbind, lapply(rows, function(r) {
    r <- r[nms]; names(r) <- nms
    as.data.frame(lapply(r, function(x) if (is.null(x)) NA else x),
                  stringsAsFactors = FALSE)
  }))
  chr <- c("arm", "bound", "status", "unbounded", "vc_names", "vc_sd",
           "msg", "file", "pdHess", "fit_sv", "mtime")
  for (nm in setdiff(nms, chr)) {
    d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
  }
  d
}

wil <- function(k, n, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  c(((p + z^2 / (2 * n)) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n),
    ((p + z^2 / (2 * n)) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n))
}

# coverage, se/sd and bias/sd for one coefficient, on a subset
cov1 <- function(tag, est, lo, hi, truth, idx = NULL) {
  if (!is.null(idx)) { est <- est[idx]; lo <- lo[idx]; hi <- hi[idx] }
  ok <- is.finite(lo) & is.finite(hi)
  est <- est[ok]; lo <- lo[ok]; hi <- hi[ok]
  n <- length(est)
  k <- sum(lo <= truth & hi >= truth)
  se <- (hi - lo) / (2 * stats::qnorm(0.975))
  w <- wil(k, n)
  cat(sprintf(paste0("  %-16s %3d/%3d = %5.1f%% (%4.1f,%5.1f)  se/sd ",
                     "%6.3f  sd %8.5f  meanse %9.5f  bias/sd %6.3f\n"),
              tag, k, n, 100 * k / n, 100 * w[1], 100 * w[2],
              mean(se) / stats::sd(est), stats::sd(est), mean(se),
              (mean(est) - truth) / stats::sd(est)))
  invisible(c(k = k, n = n, sesd = mean(se) / stats::sd(est)))
}

root <- "dev/eamhier-rec"
A <- rev_read(file.path(root, "A"))
C <- rev_read(file.path(root, "C"))
S <- rev_read(file.path(root, "S1"))
B <- rev_read(file.path(root, "B"))
cat("records: A", nrow(A), " C", nrow(C), " S1", nrow(S),
    " B", nrow(B), "\n")
cat("status ok: A", sum(A$status == "ok"), " C", sum(C$status == "ok"),
    " S1", sum(S$status == "ok"), " B", sum(B$status == "ok"), "\n\n")

tr <- list(mu0 = 0.4, mu_cond = 0.9, bs = 1.4, ndt = 0.25,
           ndt_mean = 0.25 * exp(0.12^2 / 2), sd_mu = 0.35,
           sd_lbs = 0.20, sv = 0.4)

## ---------------------------------------------------------------- 2.
## THE RETRACTION. Arm C at 60 and at the first 14, same statistic.
cat("== ARM C, all 60 ==\n")
cov1("mu0", C$mu0, C$mu0_lo, C$mu0_hi, tr$mu0)
cov1("mu_cond", C$mu_cond, C$mu_cond_lo, C$mu_cond_hi, tr$mu_cond)
cov1("log bs", C$lbs, C$lbs_lo, C$lbs_hi, log(tr$bs))
cov1("log sv", C$lsv, C$lsv_lo, C$lsv_hi, log(tr$sv))
cov1("sd(mu|s)", C$sd1, C$sd1_lo, C$sd1_hi, tr$sd_mu)
cov1("sd(log bs|s)", C$sd2, C$sd2_lo, C$sd2_hi, tr$sd_lbs)
cov1("sd(ndt|s)", C$sd3, C$sd3_lo, C$sd3_hi, NA_real_)

ord_seed <- order(C$seed)
ord_time <- order(C$mtime)
for (nm in c("seed", "mtime")) {
  idx <- if (nm == "seed") ord_seed[1:14] else ord_time[1:14]
  cat("\n== ARM C, FIRST 14 by ", nm, ": seeds ",
      paste(range(C$seed[idx]), collapse = " to "), " ==\n", sep = "")
  cov1("mu_cond", C$mu_cond, C$mu_cond_lo, C$mu_cond_hi, tr$mu_cond,
       idx)
  cov1("log sv", C$lsv, C$lsv_lo, C$lsv_hi, log(tr$sv), idx)
}

cat("\n-- se/sd for log sv over a growing prefix of arm C, by seed --\n")
for (n in c(10, 12, 13, 14, 15, 16, 18, 20, 25, 30, 40, 50, 60)) {
  idx <- ord_seed[seq_len(n)]
  e <- C$lsv[idx]; lo <- C$lsv_lo[idx]; hi <- C$lsv_hi[idx]
  se <- (hi - lo) / (2 * stats::qnorm(0.975))
  e2 <- C$mu_cond[idx]
  se2 <- (C$mu_cond_hi[idx] - C$mu_cond_lo[idx]) /
    (2 * stats::qnorm(0.975))
  cat(sprintf("  n=%2d  log sv se/sd %6.3f   mu_cond se/sd %6.3f\n",
              n, mean(se) / stats::sd(e), mean(se2) / stats::sd(e2)))
}

cat("\n-- the same growing prefix, arm A (no sv) --\n")
oa <- order(A$seed)
for (n in c(14, 30, 60)) {
  idx <- oa[seq_len(n)]
  se2 <- (A$mu_cond_hi[idx] - A$mu_cond_lo[idx]) /
    (2 * stats::qnorm(0.975))
  cat(sprintf("  n=%2d  mu_cond se/sd %6.3f\n", n,
              mean(se2) / stats::sd(A$mu_cond[idx])))
}

cat("\n== ARM A, all 60 ==\n")
cov1("mu0", A$mu0, A$mu0_lo, A$mu0_hi, tr$mu0)
cov1("mu_cond", A$mu_cond, A$mu_cond_lo, A$mu_cond_hi, tr$mu_cond)
cov1("log bs", A$lbs, A$lbs_lo, A$lbs_hi, log(tr$bs))
cov1("sd(mu|s)", A$sd1, A$sd1_lo, A$sd1_hi, tr$sd_mu)
cov1("sd(log bs|s)", A$sd2, A$sd2_lo, A$sd2_hi, tr$sd_lbs)
cat(sprintf("  sd(log bs|s) mean %.6f mcse %.6f z %.3f n %d\n",
            mean(A$sd2), stats::sd(A$sd2) / sqrt(nrow(A)),
            (mean(A$sd2) - 0.2) / (stats::sd(A$sd2) / sqrt(nrow(A))),
            nrow(A)))
cat(sprintf("  sd(mu|s)     mean %.6f mcse %.6f z %.3f\n",
            mean(A$sd1), stats::sd(A$sd1) / sqrt(nrow(A)),
            (mean(A$sd1) - 0.35) / (stats::sd(A$sd1) / sqrt(nrow(A)))))
cat(sprintf("  bs           mean %.6f mcse %.6f z %.3f\n",
            mean(exp(A$lbs)), stats::sd(exp(A$lbs)) / sqrt(nrow(A)),
            (mean(exp(A$lbs)) - 1.4) /
              (stats::sd(exp(A$lbs)) / sqrt(nrow(A)))))
cat(sprintf("  ndt rmse ms  mean %.4f   bias ms mean %.4f\n",
            mean(A$ndt_rmse_ms), mean(A$ndt_bias_ms)))
cat(sprintf("  diagnose unbounded_dpar values: %s\n",
            paste(names(table(A$unbounded)), table(A$unbounded),
                  sep = ":", collapse = " ")))
cat(sprintf("  C unbounded_dpar values: %s\n",
            paste(names(table(C$unbounded)), table(C$unbounded),
                  sep = ":", collapse = " ")))

## ---------------------------------------------------------------- 5.
## THE DRAW. ndt_pop_true recomputed from the seed with base R only.
cat("\n== 5. the draw's own mean, recomputed from the seed ==\n")
draw_mean <- function(seed, ns = 30L) {
  set.seed(seed)
  stats::rnorm(ns, 0, 0.35)          # u_mu
  stats::rnorm(ns, 0, 0.20)          # u_bs
  u <- stats::rnorm(ns, 0, 0.12)     # u_ndt
  c(mean = mean(0.25 * exp(u)), sd = stats::sd(0.25 * exp(u)))
}
dm <- t(vapply(A$seed, draw_mean, numeric(2)))
cat(sprintf("  max |record ndt_pop_true - recomputed| = %.3e\n",
            max(abs(A$ndt_pop_true - dm[, "mean"]))))
cat(sprintf("  same for arm C: %.3e\n",
            max(abs(C$ndt_pop_true -
                      t(vapply(C$seed, draw_mean,
                               numeric(2)))[, "mean"]))))
cat(sprintf("  population mean 0.25*exp(0.12^2/2) = %.6f, median 0.25\n",
            0.25 * exp(0.12^2 / 2)))
cat(sprintf("  theoretical sd of one subject's ndt = %.6f s\n",
            0.25 * exp(0.0072) * sqrt(exp(0.0144) - 1)))
cat(sprintf("  theoretical se of a 30-subject mean  = %.6f s (%.2f ms)\n",
            0.25 * exp(0.0072) * sqrt(exp(0.0144) - 1) / sqrt(30),
            1000 * 0.25 * exp(0.0072) * sqrt(exp(0.0144) - 1) /
              sqrt(30)))
cat(sprintf("  observed sd of the 60 drawn means    = %.6f s (%.2f ms)\n",
            stats::sd(dm[, "mean"]), 1000 * stats::sd(dm[, "mean"])))
cat(sprintf("  mean reported ndt_se over arm A      = %.6f s (%.2f ms)\n",
            mean(A$ndt_se), 1000 * mean(A$ndt_se)))
zfix <- abs(A$ndt_hat - 0.25) / A$ndt_se
zdrw <- abs(A$ndt_hat - A$ndt_pop_true) / A$ndt_se
cat(sprintf(paste0("  arm A: z<4 against 0.25 on %d/60, max ",
                   "%.2f | against the draw %d/60, max %.2f\n"),
            sum(zfix < 4), max(zfix), sum(zdrw < 4), max(zdrw)))
zfixC <- abs(C$ndt_hat - 0.25) / C$ndt_se
zdrwC <- abs(C$ndt_hat - C$ndt_pop_true) / C$ndt_se
cat(sprintf(paste0("  arm C: z<4 against 0.25 on %d/60, max ",
                   "%.2f | against the draw %d/60, max %.2f\n"),
            sum(zfixC < 4), max(zfixC), sum(zdrwC < 4), max(zdrwC)))
cat("  arm A seeds failing z<4 against 0.25:\n")
j <- which(zfix >= 4)
for (i in j) {
  cat(sprintf("    seed %d  z %6.2f  ndt_hat %.6f  draw mean %.6f  se %.6f\n",
              A$seed[i], zfix[i], A$ndt_hat[i], A$ndt_pop_true[i],
              A$ndt_se[i]))
}
# the predicted failure rate of the assertion, from arithmetic alone
sd_draw <- stats::sd(dm[, "mean"])
sd_tot <- sqrt(sd_draw^2 + mean(A$ndt_se)^2)
cat(sprintf(paste0("  predicted P(|z|>4) = ",
                   "2*pnorm(-4*%.6f/%.6f) = %.4f -> %.1f of 60\n"),
            mean(A$ndt_se), sd_tot,
            2 * stats::pnorm(-4 * mean(A$ndt_se) / sd_tot),
            60 * 2 * stats::pnorm(-4 * mean(A$ndt_se) / sd_tot)))
cat(sprintf(paste0("  rmse/spread under 1 on %d of 60 (A, ",
                   "max %.3f), %d of 60 (C, max %.3f)\n"),
            sum(A$ndt_rmse_ms / A$ndt_sd_true_ms < 1),
            max(A$ndt_rmse_ms / A$ndt_sd_true_ms),
            sum(C$ndt_rmse_ms / C$ndt_sd_true_ms < 1),
            max(C$ndt_rmse_ms / C$ndt_sd_true_ms)))

## ---------------------------------------------------------------- 3.
## ARM B, the "1.0000 sd 0.0000" ratio.
cat("\n== 3. arm B: sd(ndt) / (fraction x sd(floor)) ==\n")
pg <- B[B$bound == "pg" & B$status == "ok", , drop = FALSE]
gl <- B[B$bound == "gl" & B$status == "ok", , drop = FALSE]
j <- match(pg$seed, gl$seed)
cat("  paired seeds:", sum(!is.na(j)), "\n")
cat("  pg vc_names:", pg$vc_names[1], "\n")
r <- pg$ndt_sd_hat_ms / (pg$ndt_frac * pg$floor_sd_ms)
cat(sprintf("  ratio: mean %.10f  sd %.3e  min %.10f  max %.10f\n",
            mean(r), stats::sd(r), min(r), max(r)))
cat(sprintf("  max |ratio - 1| = %.3e ; ratio - 1 in ulps of 1: %.1f\n",
            max(abs(r - 1)), max(abs(r - 1)) / .Machine$double.eps))
cat(sprintf("  pg ndt variance component sd2: mean %.3e max %.3e\n",
            mean(pg$sd2), max(pg$sd2)))
cat(sprintf("  pg sd2 < 1e-4 on %d of %d\n", sum(pg$sd2 < 1e-4),
            nrow(pg)))
cat(sprintf("  correlation of ratio-1 with sd2: %.4f\n",
            suppressWarnings(stats::cor(r - 1, pg$sd2))))
cat(sprintf("  glb arm: sd(ndt) mean %.4f ms, its ratio mean %.6f sd %.3e\n",
            mean(gl$ndt_sd_hat_ms),
            mean(gl$ndt_sd_hat_ms / (gl$ndt_frac * gl$floor_sd_ms)),
            stats::sd(gl$ndt_sd_hat_ms /
                        (gl$ndt_frac * gl$floor_sd_ms))))
cat(sprintf("  paired ndt bias ms: pg %.4f gl %.4f diff %.4f +- %.4f\n",
            mean(1000 * (pg$ndt_hat - 0.25)),
            mean(1000 * (gl$ndt_hat[j] - 0.25)),
            mean(1000 * (pg$ndt_hat - gl$ndt_hat[j])),
            stats::sd(1000 * (pg$ndt_hat - gl$ndt_hat[j])) /
              sqrt(nrow(pg))))
cat(sprintf("  paired sd(ndt) ms : pg %.4f gl %.4f diff %.4f +- %.4f\n",
            mean(pg$ndt_sd_hat_ms), mean(gl$ndt_sd_hat_ms[j]),
            mean(pg$ndt_sd_hat_ms - gl$ndt_sd_hat_ms[j]),
            stats::sd(pg$ndt_sd_hat_ms - gl$ndt_sd_hat_ms[j]) /
              sqrt(nrow(pg))))

## ---------------------------------------------------------------- 1.
## S1 and the collapse.
cat("\n== 1. arm S1 ==\n")
cat(sprintf("  fit_sv TRUE on %d, sv_true 0.4 on %d\n",
            sum(S$fit_sv == "TRUE"), sum(S$sv_true == 0.4)))
sv <- exp(S$lsv)
o <- order(sv)
cat(sprintf("  sv range %.6g to %.6g; smallest three: %s\n",
            min(sv), max(sv),
            paste(sprintf("%.4g (seed %d)", sv[o][1:3], S$seed[o][1:3]),
                  collapse = ", ")))
se_lsv <- (S$lsv_hi - S$lsv_lo) / (2 * stats::qnorm(0.975))
cat(sprintf("  reported se(log sv): median %.5f max %.2f ratio %.1f\n",
            stats::median(se_lsv), max(se_lsv),
            max(se_lsv) / stats::median(se_lsv)))
k <- which.max(se_lsv)
cat(sprintf(paste0("  the outlier: seed %d  log sv %.4f  ci (%.2f, ",
                   "%.2f)  conv %d  pdHess %s  n_bad_se %d\n"),
            S$seed[k], S$lsv[k], S$lsv_lo[k], S$lsv_hi[k], S$conv[k],
            S$pdHess[k], S$n_bad_se[k]))
cov1("log sv all 60", S$lsv, S$lsv_lo, S$lsv_hi, log(0.4))
cov1("log sv drop 1", S$lsv[-k], S$lsv_lo[-k], S$lsv_hi[-k], log(0.4))
cov1("mu_cond", S$mu_cond, S$mu_cond_lo, S$mu_cond_hi, 0.9)
miss <- which(!(S$mu_cond_lo <= 0.9 & S$mu_cond_hi >= 0.9))
cat(sprintf(paste0("  mu_cond misses at seeds %s, their sv ",
                   "%s (mean %.4f); covering mean %.4f\n"),
            paste(S$seed[miss], collapse = " "),
            paste(sprintf("%.3f", sv[miss]), collapse = " "),
            mean(sv[miss]), mean(sv[-miss])))
cat(sprintf(paste0("  Spearman(log sv, mu_cond) = %.4f ",
                   "(all 60), %.4f (drop the collapse)\n"),
            stats::cor(S$lsv, S$mu_cond, method = "spearman"),
            stats::cor(S$lsv[-k], S$mu_cond[-k], method = "spearman")))

## ---------------------------------------------------------------- 6.
## POWER.
cat("\n== 6. power at n = 60 ==\n")
excl <- vapply(0:60, function(k) {
  w <- wil(k, 60); w[2] < 0.95 || w[1] > 0.95
}, logical(1))
cat("  counts whose Wilson interval excludes 0.95:",
    paste(range((0:60)[excl & (0:60) > 30]), collapse = " to "), "\n")
cat("  largest such count:", max((0:60)[excl & (0:60) < 58]), "\n")
pw <- stats::pbinom(53, 60, 0.90)
cat(sprintf("  P(X <= 53 | n=60, p=0.90) = %.4f  (Wilson-rule power)\n",
            pw))
cat(sprintf("  P(X <= 53 | n=60, p=0.95) = %.4f  (its size)\n",
            stats::pbinom(53, 60, 0.95)))
# exact binomial test rejection region at alpha 0.05, p0 = 0.95
rej <- vapply(0:60, function(k)
  stats::binom.test(k, 60, 0.95)$p.value < 0.05, logical(1))
kk <- max((0:60)[rej & (0:60) < 57])
cat(sprintf(paste0("  exact binom.test rejects at X <= %d; ",
                   "power at p=0.90 = %.4f, size %.4f\n"),
            kk, stats::pbinom(kk, 60, 0.90),
            stats::pbinom(kk, 60, 0.95)))
# normal-approximation one-sample proportion power, two-sided
za <- stats::qnorm(0.975)
p0 <- 0.95; p1 <- 0.90
pw_norm <- stats::pnorm((abs(p1 - p0) * sqrt(60) -
                           za * sqrt(p0 * (1 - p0))) /
                          sqrt(p1 * (1 - p1)))
cat(sprintf("  normal-approx power (se under H0) = %.4f\n", pw_norm))
pw_norm2 <- stats::pnorm((abs(p1 - p0) * sqrt(60) -
                            za * sqrt(p1 * (1 - p1))) /
                           sqrt(p1 * (1 - p1)))
cat(sprintf("  normal-approx power (se under H1) = %.4f\n", pw_norm2))
n_need <- (za * sqrt(p0 * (1 - p0)) +
             stats::qnorm(0.8) * sqrt(p1 * (1 - p1)))^2 / (p1 - p0)^2
cat(sprintf("  n for 80%% power, normal approx = %.1f\n", n_need))
for (n in c(180, 200, 250, 300, 400, 435, 500)) {
  kk <- max((0:n)[vapply(0:n, function(k)
    stats::binom.test(k, n, 0.95)$p.value < 0.05,
    logical(1)) & (0:n) < 0.95 * n])
  cat(sprintf(paste0("  n=%3d: exact test rejects at X<=%d, power at ",
                     "p=0.90 = %.4f | Wilson-rule power %.4f\n"),
              n, kk, stats::pbinom(kk, n, 0.90),
              stats::pbinom(max((0:n)[vapply(0:n, function(k) {
                w <- wil(k, n); w[2] < 0.95
              }, logical(1))]), n, 0.90)))
}
cat(sprintf("  sd of a coverage near 0.95 at n=60: %.4f (%.2f points)\n",
            sqrt(0.95 * 0.05 / 60), 100 * sqrt(0.95 * 0.05 / 60)))

## ------------------------------------------------------------- cost.
cat("\n== cost ==\n")
for (nm in c("A", "C")) {
  x <- if (nm == "A") A else C
  cat(sprintf(paste0("  arm %s fit_s: median %.1f max %.1f min ",
                     "%.1f | peak_ws median %.0f max %.0f MB\n"),
              nm, stats::median(x$fit_s), max(x$fit_s), min(x$fit_s),
              stats::median(x$peak_ws_mb), max(x$peak_ws_mb)))
}
cat(sprintf("  S1 fit_s: median %.1f max %.1f\n",
            stats::median(S$fit_s), max(S$fit_s)))
cat(sprintf("  B  fit_s: median %.1f max %.1f\n",
            stats::median(B$fit_s), max(B$fit_s)))
cat(sprintf("  record mtime span: A %s to %s\n", min(A$mtime),
            max(A$mtime)))
cat(sprintf("                     C %s to %s\n", min(C$mtime),
            max(C$mtime)))
cat(sprintf("                     S1 %s to %s\n", min(S$mtime),
            max(S$mtime)))
cat(sprintf("                     B %s to %s\n", min(B$mtime),
            max(B$mtime)))
