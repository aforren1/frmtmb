# Lane `latent`, item 2.4 of dev/extension-gaps-plan.md.
#
# lca() at the plan's realistic scale, K = 4 and n = 2000 with ten
# binary items and two covariates on membership, scored against the
# simulator's truth and against poLCA on the SAME data.
#
# Two questions, not one:
#   recovery  does the family get the truth back, and does its own
#             interval cover it at the stated rate;
#   identity  does poLCA's EM reach the same optimum from its own
#             starts as lca()'s single deterministic start does.
#
# The plan says "none expected" for this row. That is what is tested
# here rather than assumed.
#
#   Rscript dev/latent-2p4-lca.R [reps] [out.tsv]
#
# Seeds are 20260909 + r for replicate r; the poLCA arm's own start
# seed is 700000 + r. Both are written into every output row.

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("dev/latent-lca-sim.R")
suppressMessages(loadNamespace("poLCA"))

args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) >= 1L) as.integer(args[[1L]]) else 100L
OUT <- if (length(args) >= 2L) args[[2L]] else "dev/latent-2p4-lca.tsv"
NREP_POLCA <- 10L
K <- 4L

# poLCA agreement needs the optimizer run to the depth its EM reaches
# at tol = 1e-12; this is test-lca.R's own `tight()`.
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))

elapsed <- function(expr) {
  t0 <- Sys.time()
  force(expr)
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

# The K x J endorsement table an lca() fit implies: probability of
# category 2 (the "yes" code the simulator draws) for each class and
# item.
frm_profile_matrix <- function(fit) {
  pf <- lca_profiles(fit)
  vapply(pf, function(m) m[, 2L], numeric(nrow(pf[[1L]])))
}

polca_profile_matrix <- function(pl) {
  vapply(pl$probs, function(m) m[, 2L], numeric(nrow(pl$probs[[1L]])))
}

# The gating coefficients as a K x p matrix on a NAMED reference class.
# frmtmb's theta dpars are class k against class K, so row K is zero;
# `ref` re-references to whichever row the comparison wants.
frm_gating_matrix <- function(fit, K) {
  b <- fixef(fit)
  p <- length(b[[1L]])
  B <- matrix(0, K, p, dimnames = list(NULL, names(b[[1L]])))
  for (k in seq_len(K - 1L)) B[k, ] <- unname(b[[paste0("theta", k)]])
  B
}

# Var(theta_a[j] - theta_b[j]) from the fit's own covariance, with a
# zero row for class K, which is fixed rather than estimated.
frm_gating_contrast_se <- function(fit, a, b, j, K) {
  V <- vcov(fit)
  nms <- rownames(V)
  p <- length(fixef(fit)[[1L]])
  cn <- names(fixef(fit)[[1L]])[j]
  L <- rep(0, length(nms))
  hit <- function(k, sgn) {
    if (k >= K) return(invisible(NULL))
    i <- which(nms == paste0("theta", k, "_", cn))
    if (length(i) != 1L) stop("no vcov row for theta", k, "_", cn)
    L[i] <<- L[i] + sgn
  }
  hit(a, 1); hit(b, -1)
  sqrt(drop(crossprod(L, V %*% L)))
}

cf_names <- as.vector(t(outer(2:4, 1:3,
                              function(i, j) paste0("c", i, "_", j))))
hdr <- c("rep", "seed", "seed_polca", "perm_truth", "ll_frm", "ll_polca",
         "ll_rel_gap", "frm_beats_polca", "prof_id_max",
         "prof_id_over_polca_se", "post_id_max", "coef_id_max",
         "coef_id_over_polca_se", "prof_err_max", "prof_err_rmse",
         "coef_err_max", "coef_cover_9", "share_err_max", "entropy",
         "modal_acc", "s_frm", "s_polca", "maxgrad_rel", "pdhess",
         paste0("err_", cf_names), paste0("se_", cf_names))
if (!file.exists(OUT)) {
  cat(paste(hdr, collapse = "\t"), "\n", sep = "", file = OUT)
}

for (r in seq_len(REPS)) {
  seed <- 20260909L + r
  seed_pl <- 700000L + r
  s <- lca_sim(seed = seed)

  fit <- NULL
  s_frm <- elapsed(fit <- suppressWarnings(suppressMessages(
    frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
        control = tight))))
  ll_frm <- as.numeric(logLik(fit))

  dp <- data.frame(x1 = s$dd$x1, x2 = s$dd$x2)
  for (j in seq_len(ncol(s$Y))) dp[[paste0("I", j)]] <- s$Y[, j]
  ff <- stats::as.formula(paste0("cbind(",
    paste(paste0("I", seq_len(ncol(s$Y))), collapse = ", "),
    ") ~ x1 + x2"))
  set.seed(seed_pl)
  pl <- NULL
  s_pl <- elapsed(pl <- poLCA::poLCA(ff, dp, nclass = K,
                                     nrep = NREP_POLCA, verbose = FALSE,
                                     maxiter = 20000, tol = 1e-12))

  Pf <- frm_profile_matrix(fit)
  Pp <- polca_profile_matrix(pl)

  ## ---- identity against poLCA, on the same data ---------------------
  al <- lca_align(Pf, Pp)                # frm class al$perm[i] ~ pl i
  prof_id <- al$dist
  pl_se <- vapply(pl$probs.se, function(m) m[, 2L],
                  numeric(nrow(pl$probs.se[[1L]])))
  prof_id_rel <- prof_id / max(pl_se)

  post_f <- lca_probs(fit)
  post_id <- max(abs(post_f[, al$perm, drop = FALSE] - pl$posterior))

  Bf <- frm_gating_matrix(fit, K)[al$perm, , drop = FALSE]
  Bf <- sweep(Bf, 2L, Bf[1L, ], "-")
  Bp <- rbind(0, t(pl$coeff))
  coef_id <- max(abs(Bf - Bp))
  pl_cse <- rbind(NA_real_, t(pl$coeff.se))
  coef_id_rel <- coef_id / max(pl_cse[-1L, ])

  ## ---- recovery against the simulator's truth -----------------------
  alt <- lca_align(Pf, s$base)           # frm class alt$perm[i] ~ true i
  prof_err <- alt$dist
  prof_rmse <- sqrt(mean((Pf[alt$perm, , drop = FALSE] - s$base)^2))
  Bt <- frm_gating_matrix(fit, K)[alt$perm, , drop = FALSE]
  Bt <- sweep(Bt, 2L, Bt[1L, ], "-")
  coef_err <- max(abs(Bt - s$gam))
  cf_err <- numeric(0)
  cf_se <- numeric(0)
  for (i in 2:K) {
    for (j in seq_len(ncol(Bt))) {
      cf_err <- c(cf_err, Bt[i, j] - s$gam[i, j])
      cf_se <- c(cf_se, frm_gating_contrast_se(fit, alt$perm[i],
                                               alt$perm[1L], j, K))
    }
  }
  cover <- sum(abs(cf_err) < stats::qnorm(0.975) * cf_se)
  sizes <- attr(lca_profiles(fit), "class_sizes")[alt$perm]
  share_err <- max(abs(unname(sizes) - unname(prop.table(table(s$cl)))))
  ent <- attr(post_f, "entropy")
  acc <- mean(max.col(post_f[, alt$perm, drop = FALSE]) == s$cl)

  d1 <- frmtmb::diagnose(fit, quiet = TRUE)

  row <- c(r, seed, seed_pl, paste(alt$perm, collapse = ""),
           formatC(ll_frm, digits = 12, format = "f"),
           formatC(pl$llik, digits = 12, format = "f"),
           formatC(abs(ll_frm - pl$llik) / abs(ll_frm), digits = 4,
                   format = "e"),
           as.integer(ll_frm - pl$llik > 1e-9 * abs(ll_frm)),
           formatC(c(prof_id, prof_id_rel, post_id, coef_id, coef_id_rel,
                     prof_err, prof_rmse, coef_err), digits = 6,
                   format = "g"),
           cover,
           formatC(c(share_err, ent, acc, s_frm, s_pl,
                     d1$max_grad / abs(ll_frm)), digits = 6, format = "g"),
           isTRUE(d1$pdHess),
           formatC(c(cf_err, cf_se), digits = 6, format = "g"))
  cat(paste(row, collapse = "\t"), "\n", sep = "", file = OUT,
      append = TRUE)
  cat(sprintf("rep %3d  dll %+0.3e  profid %.2e  proferr %.4f  cover %d/9\n",
              r, ll_frm - pl$llik, prof_id, prof_err, cover))
  flush(stdout())
}
