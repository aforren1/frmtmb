# Lane `latent`, item 2.4: does the recipe `?lca` already prints work?
#
# The "Labeling and starting values" section of `?lca` tells a user to
# do what `poLCA(nrep = 10)` does, and shows this:
#
#   p0 <- fit$frame$par_template[fit$frame$extra_names]
#   refits <- replicate(10, simplify = FALSE,
#     frm(..., start = lapply(p0, function(v) v + rnorm(length(v)))))
#
# That perturbs the item-profile logits around the DEFAULT START, which
# is not what dev/latent-2p4-modes.R tried: that one jittered around the
# OPTIMUM by two standard errors, the hmm_starts() recipe, and it never
# reached poLCA's answer. Advice on a help page is a claim, so this
# runs the page's own recipe on the eight datasets where the shipped
# start goes wrong.
#
#   Rscript dev/latent-2p4-docrecipe.R [nrep] [seeds...]

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("dev/latent-lca-sim.R")
suppressMessages(loadNamespace("poLCA"))

args <- commandArgs(trailingOnly = TRUE)
NREP <- if (length(args) >= 1L) as.integer(args[[1L]]) else 10L
SEEDS <- if (length(args) > 1L) as.integer(args[-1L]) else
  c(20260970L, 20260990L, 20260996L, 20260999L,
    20261010L, 20261013L, 20261042L, 20261043L)
K <- 4L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))

out <- "dev/latent-2p4-docrecipe.tsv"
cat("seed\tnrep\tll_cold\tll_polca10\tll_doc\tgap_cold\tgap_doc",
    "\tn_better\tn_error\tspread\tseconds\n", sep = "", file = out)

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
  set.seed(700000L + (seed %% 1000L))
  pl <- poLCA::poLCA(ff, dp, nclass = K, nrep = 10L, verbose = FALSE,
                     maxiter = 20000, tol = 1e-12)

  ## the help page's own recipe, verbatim in spirit: perturb the ITEM
  ## parameters of the DEFAULT start, refit, keep the best
  p0 <- fit$frame$par_template[fit$frame$extra_names]
  set.seed(500000L + (seed %% 1000L))
  lls <- rep(NA_real_, NREP)
  nerr <- 0L
  t0 <- Sys.time()
  for (i in seq_len(NREP)) {
    st <- lapply(p0, function(v) v + stats::rnorm(length(v)))
    f2 <- try(suppressWarnings(suppressMessages(
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
          control = tight, start = st))), silent = TRUE)
    if (inherits(f2, "try-error")) {
      nerr <- nerr + 1L
    } else {
      lls[i] <- as.numeric(logLik(f2))
    }
  }
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  best <- if (all(is.na(lls))) ll else max(c(ll, lls), na.rm = TRUE)
  nb <- sum(lls > ll + 1e-9 * abs(ll), na.rm = TRUE)
  sp <- if (sum(!is.na(lls)) > 1L) diff(range(lls, na.rm = TRUE)) else 0

  cat(paste(c(seed, NREP,
              formatC(c(ll, pl$llik, best), digits = 6, format = "f"),
              formatC(c(ll - pl$llik, best - pl$llik, nb, nerr, sp,
                        secs), digits = 6, format = "g")),
            collapse = "\t"), "\n", sep = "", file = out, append = TRUE)
  cat(sprintf("seed %d: cold %.3f  poLCA10 %.3f  docrecipe %.3f",
              seed, ll, pl$llik, best),
      sprintf("  gap_cold %+.2f gap_doc %+.2f  better %d/%d err %d",
              ll - pl$llik, best - pl$llik, nb, NREP, nerr),
      sprintf("  spread %.2f  %.1fs\n", sp, secs))
  flush(stdout())
}
