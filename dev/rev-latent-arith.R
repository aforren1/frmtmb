# Reviewer: recompute every count, coverage and Monte Carlo interval in
# dev/latent-findings.md straight from the lane's own output tables.
# Nothing is fitted here. The point is arithmetic and denominators: a
# coverage rate is only a rate if nothing was silently dropped from the
# bottom of the fraction.

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
rd <- function(f) utils::read.delim(paste0(D, f), stringsAsFactors = FALSE)

hr <- function(s) cat("\n==== ", s, " ====\n", sep = "")

mci <- function(k, n) {
  p <- k / n
  se <- sqrt(p * (1 - p) / n)
  sprintf("%.4f (%d/%d)  MC 95%%: %.4f to %.4f  (+/- %.2f points)",
          p, k, n, p - 1.96 * se, p + 1.96 * se, 100 * 1.96 * se)
}

## ---------------------------------------------------------------- lca
hr("lca: replicate counts and the split")
a <- rd("latent-2p4-lca.tsv")
ag <- rd("latent-2p4-lca-agree.tsv")
dg <- rd("latent-2p4-lca-disagree.tsv")
cat("rows all/agree/disagree:", nrow(a), nrow(ag), nrow(dg), "\n")
cat("seeds all: min", min(a$seed), "max", max(a$seed),
    "distinct", length(unique(a$seed)), "\n")
cat("seeds contiguous 20260910..20261109:",
    identical(sort(a$seed), 20260910:20261109), "\n")
cat("polca seeds 700001..700200:",
    identical(sort(a$seed_polca), 700001:700200), "\n")
cat("agree + disagree == all:", nrow(ag) + nrow(dg) == nrow(a), "\n")
cat("disagree seeds:", paste(sort(dg$seed), collapse = " "), "\n")
cat("any NA in agree ll_rel_gap:", anyNA(ag$ll_rel_gap), "\n")

hr("lca: how the split was defined, checked against ll_rel_gap")
cat("agree  max |ll_rel_gap|  :", format(max(abs(ag$ll_rel_gap)), digits = 4),
    "\n")
cat("disagr min |ll_rel_gap|  :", format(min(abs(dg$ll_rel_gap)), digits = 4),
    "\n")
cat("frm_beats_polca in agree :", sum(ag$frm_beats_polca), "\n")
cat("frm_beats_polca in disagr:", sum(dg$frm_beats_polca), "\n")

hr("lca: the identity claims (worst over the 192)")
cat("ll   |d|/|ll|            :", format(max(abs(ag$ll_rel_gap)), digits = 3),
    " claim 4.27e-14\n")
cat("profiles max|d|          :", format(max(ag$prof_id_max), digits = 3),
    " claim 5.61e-07\n")
cat("posterior max|d|         :", format(max(ag$post_id_max), digits = 3),
    " claim 2.12e-06\n")
cat("gating coefs max|d|      :", format(max(ag$coef_id_max), digits = 3),
    " claim 2.22e-06\n")
cat("prof over polca se       :",
    format(max(ag$prof_id_over_polca_se), digits = 3), " claim 1.68e-05\n")
cat("coef over polca se       :",
    format(max(ag$coef_id_over_polca_se), digits = 3), " claim 1.11e-05\n")

hr("lca: coverage, recomputed")
cov9 <- ag$coef_cover_9
cat("coef_cover_9 range:", range(cov9), " NA:", sum(is.na(cov9)), "\n")
cat("overall:", mci(sum(cov9), 9 * nrow(ag)), " claim 94.6% 1634/1728\n")

hr("lca: the per-coefficient table, recomputed")
truth <- c(c2_1 = -0.4, c2_2 = 0.8, c2_3 = -0.5,
           c3_1 = 0.2, c3_2 = -0.6, c3_3 = 0.9,
           c4_1 = -0.1, c4_2 = 0.3, c4_3 = 0.4)
tab <- do.call(rbind, lapply(names(truth), function(nm) {
  e <- ag[[paste0("err_", nm)]]
  s <- ag[[paste0("se_", nm)]]
  cover <- mean(abs(e) <= stats::qnorm(0.975) * s)
  data.frame(coef = nm, truth = truth[[nm]],
             bias = mean(e), rmse = sqrt(mean(e^2)),
             empSD = stats::sd(e), meanSE = mean(s),
             cover95 = 100 * cover, nNA = sum(is.na(e) | is.na(s)))
}))
print(tab, row.names = FALSE, digits = 4)
cat("sum of per-coef covered:", sum(round(tab$cover95 / 100 * nrow(ag))),
    " vs coef_cover_9 total ", sum(cov9), "\n")

hr("lca: the other-quantity table")
qt <- function(v) sprintf("%.4f / %.4f / %.4f", min(v), stats::median(v),
                          max(v))
cat("prof_err_max  min/med/max:", qt(ag$prof_err_max),
    " claim 0.034/0.058/0.098\n")
cat("prof_err_rmse            :", qt(ag$prof_err_rmse),
    " claim 0.016/0.023/0.034\n")
cat("share_err_max            :", qt(ag$share_err_max),
    " claim 0.0027/0.011/0.031\n")
cat("entropy                  :", qt(ag$entropy),
    " claim 0.784/0.807/0.830\n")
cat("modal_acc                :", qt(ag$modal_acc),
    " claim 0.878/0.897/0.916\n")
cat("pdhess TRUE              :", sum(ag$pdhess), "of", nrow(ag), "\n")
cat("maxgrad_rel range        :", format(range(ag$maxgrad_rel), digits = 3),
    " claim 7.6e-10 to 8.0e-09\n")

hr("lca: the eight, recomputed from latent-2p4-modes.tsv")
m <- rd("latent-2p4-modes.tsv")
print(m[, c("seed", "ll_cold", "ll_polca10", "gap", "conv", "maxgrad_rel",
            "pdhess", "polca1_low", "polca1_n")], row.names = FALSE,
      digits = 10)
cat("gap range:", format(range(m$gap), digits = 6), " claim -243 to -284\n")
cat("pdHess TRUE count:", sum(m$pdhess), "of", nrow(m), "\n")
cat("maxgrad_rel range:", format(range(m$maxgrad_rel), digits = 3),
    " claim 1.9e-09 to 2.7e-08\n")
cat("nlminb code 0 count:", sum(m$conv == 0), "\n")
cat("pass pdHess AND maxgrad_rel<1e-3:",
    sum(m$pdhess & m$maxgrad_rel < 1e-3), " claim 7 of 8\n")
cat("polca1_low range:", range(m$polca1_low), "of", unique(m$polca1_n),
    " claim 0 to 4 of 10\n")
cat("jitter arm best gap range:", format(range(m$jit_gap), digits = 5),
    " claim 72 to 273 units left\n")
cat("jit_nbetter:", paste(m$jit_nbetter, collapse = " "), "\n")

hr("lca: the doc recipe")
dr <- rd("latent-2p4-docrecipe.tsv")
print(dr, row.names = FALSE, digits = 8)
cat("reached polca optimum (|gap_doc| < 1e-6):",
    sum(abs(dr$gap_doc) < 1e-6), "of", nrow(dr), "\n")
cat("gap_doc range:", format(range(abs(dr$gap_doc)), digits = 3),
    " claim 3e-11 to 6e-11\n")
cat("seconds range:", format(range(dr$seconds), digits = 4),
    " findings claim 6.1 to 16.8 s; ?lca says 'about 5 seconds'\n")
cat("n_better range:", range(dr$n_better), " claim 8 to 10 of 10\n")

## ------------------------------------------------------------- hmm A
hr("hmm arm A: counts and coverage")
A <- rd("latent-2p3-hmm-A.tsv")
cat("rows:", nrow(A), " distinct seeds:", length(unique(A$seed)),
    " contiguous 20260910..20260949:",
    identical(sort(A$seed), 20260910:20260949), "\n")
covA <- grep("^cov_", names(A), value = TRUE)
covA <- covA[covA != "cov_theta1"]
cat("coverage columns used:", length(covA), ":", paste(covA, collapse = " "),
    "\n")
cat("cov_theta1 all NA:", all(is.na(A$cov_theta1)), "\n")
mA <- as.matrix(A[, covA])
cat("NA cells:", sum(is.na(mA)), " non-0/1 cells:",
    sum(!is.na(mA) & !(mA %in% c(0, 1))), "\n")
cat("overall:", mci(sum(mA), length(mA)), " claim 95.6% 459/480\n")
truthA <- c(mu1 = -2, mu2 = 0, mu3 = 3,
            sigma1 = log(0.7), sigma2 = log(0.7), sigma3 = log(0.7),
            tr12 = -1.5, tr13 = -2.5, tr22 = 2, tr23 = -1,
            tr32 = -1, tr33 = 2)
tA <- do.call(rbind, lapply(names(truthA), function(nm) {
  e <- A[[paste0("err_", nm)]]; s <- A[[paste0("se_", nm)]]
  data.frame(par = nm, truth = truthA[[nm]], bias = mean(e),
             rmse = sqrt(mean(e^2)), empSD = stats::sd(e),
             meanSE = mean(s), cover95 = 100 * mean(A[[paste0("cov_", nm)]]))
}))
print(tA, row.names = FALSE, digits = 4)
cat("perm all identity:", paste(unique(A$perm), collapse = " | "), "\n")
cat("pdhess:", sum(A$pdhess), "of", nrow(A), "\n")
cat("maxgrad_rel range:", format(range(A$maxgrad_rel), digits = 3),
    " claim 4.3e-08 to 4.9e-07\n")
cat("max nlminb_gr:", max(A$nlminb_gr), " claim <= 42\n")
cat("mu_err_max max:", format(max(A$mu_err_max), digits = 3), " claim 0.021\n")
cat("sd_err_max max:", format(max(A$sd_err_max), digits = 3), " claim 0.017\n")
cat("tpm_err_max max:", format(max(A$tpm_err_max), digits = 3),
    " claim 0.014\n")

hr("hmm arm A: the identity, 3 replicates")
Ai <- rd("latent-2p3-hmm-A-id.tsv")
print(Ai[, c("rep", "seed", "ll_frm", "ll_ref", "ll_rel_gap", "est_max_gap",
             "tpm_max_gap", "identity_ok")], row.names = FALSE, digits = 12)
cat("ll_rel_gap range:", format(range(abs(Ai$ll_rel_gap)), digits = 3),
    " claim 3.1e-14 to 4.9e-12\n")
cat("est_max_gap range:", format(range(Ai$est_max_gap), digits = 3),
    " claim 7.9e-07 to 1.5e-06\n")
cat("tpm_max_gap range:", format(range(Ai$tpm_max_gap), digits = 3),
    " claim 4.0e-07 to 5.9e-07\n")
cat("mean se on a state mean, over the 3:",
    format(range(c(Ai$se_mu1, Ai$se_mu2, Ai$se_mu3)), digits = 3),
    " claim 0.0081 to 0.0095\n")

## ------------------------------------------------------------- hmm B
hr("hmm arm B: counts and coverage")
B <- rd("latent-2p3-hmm-B.tsv")
cat("rows:", nrow(B), " contiguous 20260910..20260924:",
    identical(sort(B$seed), 20260910:20260924), "\n")
covB <- grep("^cov_", names(B), value = TRUE)
mB <- as.matrix(B[, covB])
cat("coverage columns used:", length(covB), "\n")
cat("NA cells:", sum(is.na(mB)), " non-0/1:",
    sum(!is.na(mB) & !(mB %in% c(0, 1))), "\n")
cat("overall:", mci(sum(mB, na.rm = TRUE), sum(!is.na(mB))),
    " claim 94.4% 184/195\n")
truthB <- c(truthA, theta1 = log(0.6))
tB <- do.call(rbind, lapply(names(truthB), function(nm) {
  e <- B[[paste0("err_", nm)]]; s <- B[[paste0("se_", nm)]]
  data.frame(par = nm, truth = truthB[[nm]], bias = mean(e),
             rmse = sqrt(mean(e^2)), empSD = stats::sd(e),
             meanSE = mean(s), cover95 = 100 * mean(B[[paste0("cov_", nm)]]))
}))
print(tB, row.names = FALSE, digits = 4)
cat("sd_tr12 mean:", format(mean(B$sd_tr12), digits = 5),
    " sd:", format(stats::sd(B$sd_tr12), digits = 4),
    " range:", format(range(B$sd_tr12), digits = 4),
    " claim 0.615 +/- 0.067, 0.511 to 0.718\n")
cat("pdhess:", sum(B$pdhess), "of", nrow(B), "\n")
cat("maxgrad_rel range:", format(range(B$maxgrad_rel), digits = 3),
    " claim 4.5e-08 to 3.5e-07\n")
cat("max nlminb_gr:", max(B$nlminb_gr), " claim <= 31\n")
cat("mu_err_max max:", format(max(B$mu_err_max), digits = 3), " claim 0.017\n")
cat("sd_err_max max:", format(max(B$sd_err_max), digits = 3), " claim 0.018\n")
cat("tpm_err_max max:", format(max(B$tpm_err_max), digits = 3),
    " claim 0.030\n")
cat("identity rows (ll_ref not NA):", sum(!is.na(B$ll_ref)), " claim 6\n")
idx <- !is.na(B$ll_ref)
cat("ll_rel_gap over those:", format(range(abs(B$ll_rel_gap[idx])),
                                     digits = 3),
    " claim 6.1e-12 to 6.9e-11\n")
cat("est_max_gap over those:", format(range(B$est_max_gap[idx]), digits = 3),
    " claim 5.2e-06 to 1.3e-04\n")
cat("tpm_max_gap over those:", format(range(B$tpm_max_gap[idx]), digits = 3),
    " claim 0.039 to 0.180\n")

## ------------------------------------------------------- starts probe
hr("hmm_starts probe: the jitter sweep")
sp <- rd("latent-2p3-starts-probe.tsv")
print(sp, row.names = FALSE, digits = 10)
agg <- aggregate(cbind(recovered, converged, not_converged, error) ~ jitter,
                 data = sp, FUN = sum)
agg$reps <- as.numeric(table(sp$jitter))
print(agg, row.names = FALSE)
cat("n per row:", paste(unique(sp$n), collapse = " "), "\n")
