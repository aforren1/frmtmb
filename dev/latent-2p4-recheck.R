# Lane `latent`, punch round 1: the lca recovery table, re-measured
# against the SHIPPED default start after it changed.
#
# The 200-replicate sweep in `dev/latent-2p4-lca.R` measured the old
# score-cut start. `lca_init_extras()` now slices on the response
# pattern and shrinks the class profiles halfway to the pooled ones, so
# the recovery table on `?lca` would otherwise describe a start that no
# longer ships.
#
# THE poLCA ARM IS NOT RE-RUN. Its log-likelihood for every seed is READ
# from `dev/latent-2p4-lca.tsv`, so what this costs is 200 frmtmb fits
# and 200 `confint()` calls, about 0.5 s each, and no EM at all. Nothing
# recomputed here can move the reference.
#
#   Rscript dev/latent-2p4-recheck.R [nseed] [reference tsv] [out tsv]
#
# Default seeds are the tuning block's own, 20260910 + 0..199. Point it
# at `dev/latent-2p4-oos.tsv` to do the same on the 200 fresh seeds the
# shrink weight was SCORED on rather than tuned on; that file's poLCA
# column is read the same way.

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("dev/latent-lca-sim.R")

args <- commandArgs(trailingOnly = TRUE)
NSEED <- if (length(args) >= 1L) as.integer(args[[1L]]) else 200L
REF <- if (length(args) >= 2L) args[[2L]] else "dev/latent-2p4-lca.tsv"
OUTF <- if (length(args) >= 3L) args[[3L]] else
  "dev/latent-2p4-recheck.tsv"
K <- 4L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))

frm_profile_matrix <- function(fit) {
  pf <- lca_profiles(fit)
  vapply(pf, function(m) m[, 2L], numeric(nrow(pf[[1L]])))
}
frm_gating_matrix <- function(fit, K) {
  b <- fixef(fit)
  B <- matrix(0, K, length(b[[1L]]),
              dimnames = list(NULL, names(b[[1L]])))
  for (k in seq_len(K - 1L)) B[k, ] <- unname(b[[paste0("theta", k)]])
  B
}
frm_gating_contrast_se <- function(fit, a, b, j, K) {
  V <- vcov(fit)
  nms <- rownames(V)
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

ref <- utils::read.delim(REF, check.names = FALSE)
cf_names <- as.vector(t(outer(2:4, 1:3,
                              function(i, j) paste0("c", i, "_", j))))
out <- OUTF
cat(paste(c("rep", "seed", "perm_truth", "ll_frm", "ll_polca",
            "ll_rel_gap", "reached_polca", "prof_err_max",
            "prof_err_rmse", "coef_err_max", "coef_cover_9",
            "share_err_max", "entropy", "modal_acc", "s_frm",
            "maxgrad_rel", "pdhess",
            paste0("err_", cf_names), paste0("se_", cf_names)),
          collapse = "\t"), "\n", sep = "", file = out)

for (i in seq_len(min(NSEED, nrow(ref)))) {
  seed <- ref$seed[i]
  s <- lca_sim(seed = seed)
  t0 <- Sys.time()
  fit <- suppressWarnings(suppressMessages(
    frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
        control = tight)))
  s_frm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ll <- as.numeric(logLik(fit))
  Pf <- frm_profile_matrix(fit)
  alt <- lca_align(Pf, s$base)
  prof_err <- alt$dist
  prof_rmse <- sqrt(mean((Pf[alt$perm, , drop = FALSE] - s$base)^2))
  Bt <- frm_gating_matrix(fit, K)[alt$perm, , drop = FALSE]
  Bt <- sweep(Bt, 2L, Bt[1L, ], "-")
  cf_err <- numeric(0)
  cf_se <- numeric(0)
  for (a in 2:K) {
    for (j in seq_len(ncol(Bt))) {
      cf_err <- c(cf_err, Bt[a, j] - s$gam[a, j])
      cf_se <- c(cf_se, frm_gating_contrast_se(fit, alt$perm[a],
                                               alt$perm[1L], j, K))
    }
  }
  cover <- sum(abs(cf_err) < stats::qnorm(0.975) * cf_se)
  post_f <- lca_probs(fit)
  sizes <- attr(lca_profiles(fit), "class_sizes")[alt$perm]
  share_err <- max(abs(unname(sizes) -
                         unname(prop.table(table(s$cl)))))
  acc <- mean(max.col(post_f[, alt$perm, drop = FALSE]) == s$cl)
  d1 <- frmtmb::diagnose(fit, quiet = TRUE)
  gapr <- abs(ll - ref$ll_polca[i]) / abs(ll)
  cat(paste(c(i, seed, paste(alt$perm, collapse = ""),
              formatC(c(ll, ref$ll_polca[i]), digits = 12,
                      format = "f"),
              formatC(gapr, digits = 4, format = "e"),
              as.integer(ll > ref$ll_polca[i] - 1e-6 * abs(ll)),
              formatC(c(prof_err, prof_rmse, max(abs(cf_err))),
                      digits = 6, format = "g"),
              cover,
              formatC(c(share_err, attr(post_f, "entropy"), acc, s_frm,
                        d1$max_grad / abs(ll)), digits = 6,
                      format = "g"),
              isTRUE(d1$pdHess),
              formatC(c(cf_err, cf_se), digits = 6, format = "g")),
            collapse = "\t"), "\n", sep = "", file = out, append = TRUE)
  if (i %% 25L == 0L) { cat("  ", i, " seeds\n", sep = ""); flush(stdout()) }
}
cat("wrote", out, "\n")
