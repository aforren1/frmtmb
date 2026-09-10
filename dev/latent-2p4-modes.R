# Lane `latent`, item 2.4: the eight replicates where lca() and poLCA
# disagreed, looked at one at a time.
#
# The 200-replicate run (dev/latent-2p4-lca.tsv) found lca()'s single
# deterministic start converging 242 to 284 log-likelihood units below
# poLCA's best of ten EM starts on 8 of 200 datasets. The plan's row
# 2.4 says "none expected", so this asks three questions the summary
# cannot answer:
#
#   1. Does anything in the fit SAY it is wrong? (diagnose(), the
#      gradient, the Hessian.)
#   2. Is it lca()'s starting values, or the surface? poLCA is run
#      again with nrep = 1 from ten separate seeds: if a single EM
#      start lands low too, the surface is multimodal and one start is
#      not enough for anybody.
#   3. Would a multistart fix it? Ten refits of lca() from starts
#      jittered by two standard errors, the same recipe hmm_starts()
#      uses for hmm().
#
#   Rscript dev/latent-2p4-modes.R [seeds...]
#
# Default seeds are the eight the 200-replicate run flagged.

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("dev/latent-lca-sim.R")
suppressMessages(loadNamespace("poLCA"))

args <- commandArgs(trailingOnly = TRUE)
SEEDS <- if (length(args)) as.integer(args) else
  c(20260970L, 20260990L, 20260996L, 20260999L,
    20261010L, 20261013L, 20261042L, 20261043L)
K <- 4L
NJIT <- 10L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))

out <- "dev/latent-2p4-modes.tsv"
cat("seed\tll_cold\tll_polca10\tgap\tconv\tmaxgrad_rel\tpdhess\tnbadse",
    "\tnflat\tnmsg\tpolca1_low\tpolca1_n\tjit_best\tjit_gap\tjit_nbetter\n",
    sep = "", file = out)

for (seed in SEEDS) {
  s <- lca_sim(seed = seed)
  dp <- data.frame(x1 = s$dd$x1, x2 = s$dd$x2)
  for (j in seq_len(ncol(s$Y))) dp[[paste0("I", j)]] <- s$Y[, j]
  ff <- stats::as.formula(paste0("cbind(",
    paste(paste0("I", seq_len(ncol(s$Y))), collapse = ", "),
    ") ~ x1 + x2"))

  fit <- suppressWarnings(suppressMessages(
    frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
        control = tight)))
  ll <- as.numeric(logLik(fit))
  d1 <- frmtmb::diagnose(fit, quiet = TRUE)
  msg <- utils::capture.output(frmtmb::diagnose(fit))

  set.seed(700000L + (seed %% 1000L))
  pl10 <- poLCA::poLCA(ff, dp, nclass = K, nrep = 10L, verbose = FALSE,
                       maxiter = 20000, tol = 1e-12)
  ll10 <- pl10$llik

  ## 2. single-start poLCA, ten separate seeds
  low <- 0L
  n1 <- 10L
  for (i in seq_len(n1)) {
    set.seed(800000L + i)
    p1 <- try(poLCA::poLCA(ff, dp, nclass = K, nrep = 1L, verbose = FALSE,
                           maxiter = 20000, tol = 1e-12), silent = TRUE)
    if (!inherits(p1, "try-error") &&
        p1$llik < ll10 - 1e-6 * abs(ll10)) {
      low <- low + 1L
    }
  }

  ## 3. a multistart of lca(), the hmm_starts() recipe by hand
  ci <- suppressWarnings(stats::confint(fit))
  par <- fit$opt$par
  se <- if (nrow(ci) == length(par)) {
    (ci[, "upr"] - ci[, "lwr"]) / (2 * stats::qnorm(0.975))
  } else rep(1, length(par))
  se[!is.finite(se) | se <= 0] <- stats::median(
    se[is.finite(se) & se > 0])
  set.seed(600000L + (seed %% 1000L))
  best <- ll
  nbetter <- 0L
  for (i in seq_len(NJIT)) {
    st <- fit$estimates
    v <- as.numeric(par) + stats::rnorm(length(par), 0, 2 * se)
    pn <- names(par)
    for (cp in unique(pn)) {
      idx <- which(pn == cp)
      if (length(st[[cp]]) == length(idx)) st[[cp]][] <- v[idx]
    }
    f2 <- try(suppressWarnings(suppressMessages(
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
          control = tight, start = st))), silent = TRUE)
    if (!inherits(f2, "try-error")) {
      l2 <- as.numeric(logLik(f2))
      if (is.finite(l2) && l2 > best + 1e-9 * abs(best)) {
        best <- l2
        nbetter <- nbetter + 1L
      }
    }
  }

  cat(paste(c(seed,
              formatC(c(ll, ll10), digits = 6, format = "f"),
              formatC(ll - ll10, digits = 6, format = "g"),
              d1$convergence,
              formatC(d1$max_grad / abs(ll), digits = 4, format = "e"),
              isTRUE(d1$pdHess), length(d1$bad_se), length(d1$flat),
              length(msg), low, n1,
              formatC(best, digits = 6, format = "f"),
              formatC(best - ll10, digits = 6, format = "g"),
              nbetter),
            collapse = "\t"), "\n", sep = "", file = out, append = TRUE)
  cat(sprintf("seed %d: cold %.4f  poLCA10 %.4f  gap %+.3f",
              seed, ll, ll10, ll - ll10),
      sprintf("  conv %d pdHess %s  poLCA nrep=1 low %d/%d",
              d1$convergence, isTRUE(d1$pdHess), low, n1),
      sprintf("  multistart best %.4f (%+.3f, %d better)\n",
              best, best - ll10, nbetter))
  cat("  diagnose() said:\n")
  cat(paste0("    ", msg, collapse = "\n"), "\n")
  flush(stdout())
}
