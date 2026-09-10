# Reviewer, punch round 1: the NEW lca start has two fallbacks of its
# own, and a new default start is a behaviour change on every lca()
# fit, not only on the design it was tuned for.
#
#   G1  `n < K`, which falls back to the score cut
#   G2  no item varies at all, which also falls back
#   G3  small and ordinary data, where the old rule was fine: does the
#       new one do worse anywhere? A remedy that fires on a correct
#       model is a real cost, and a start is a remedy that fires on
#       every fit.
#
#   Rscript dev/rev-latent-startguards.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env2.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
source(paste0(D, "latent-2p4-startfix-fns.R"))
hr <- function(s) cat("\n======== ", s, " ========\n", sep = "")
slice <- getFromNamespace("lca_init_slice", "frmtmb.latent")
extras <- getFromNamespace("lca_init_extras", "frmtmb.latent")

hr("G1 and G2: the new start's own fallbacks, by construction")
y2 <- matrix(c(1L, 2L, 1L, 2L, 2L, 1L), nrow = 2L)
cat("n = 2, K = 4 (n < K): slice =",
    paste(slice(y2, rep(2L, 3L), 4L), collapse = " "), "\n")
yc <- matrix(1L, nrow = 50L, ncol = 4L)
cat("every item constant, K = 3: slice runs and gives",
    length(unique(slice(yc, rep(2L, 4L), 3L))), "distinct labels\n")
e <- extras(yc, rep(2L, 4L), 3L)
cat("  lca_init_extras() on that returns", length(e),
    "item vectors, all finite:", all(is.finite(unlist(e))), "\n")
yna <- matrix(sample(1:2, 200, TRUE), nrow = 50L)
yna[cbind(1:10, 1:4)] <- NA_integer_
cat("with NAs present: slice runs, labels",
    paste(sort(unique(slice(yna, rep(2L, 4L), 3L))), collapse = " "),
    " extras all finite:",
    all(is.finite(unlist(extras(yna, rep(2L, 4L), 3L)))), "\n")

hr("G3: does the new start do worse anywhere the old one was fine?")
cat("Twelve small designs the old rule had no trouble with: K = 2 and\n")
cat("K = 3, n = 200 and 600, and a design where the classes DO differ\n")
cat("in how many items they endorse, which is the case the old score\n")
cat("was built for and where it should still work.\n\n")
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))
sim_mono <- function(seed, n, K, J = 6L) {
  # classes differ in HOW MANY items they endorse: the design the old
  # mean-score cut was designed for
  set.seed(seed)
  base <- matrix(0.15, K, J)
  for (k in seq_len(K)) {
    if (k > 1L) base[k, seq_len(round(J * (k - 1) / (K - 1)))] <- 0.85
  }
  cl <- sample.int(K, n, TRUE)
  Y <- matrix(0L, n, J)
  for (j in seq_len(J)) Y[, j] <- 1L + stats::rbinom(n, 1L, base[cl, j])
  dd <- data.frame(x1 = stats::rnorm(n)); dd$Y <- Y
  list(dd = dd, Y = Y)
}
sim_block <- function(seed, n, K, J = 6L) {
  set.seed(seed)
  base <- matrix(0.2, K, J)
  w <- max(1L, J %/% K)
  for (k in seq_len(K)) {
    ix <- ((k - 1L) * w + 1L):min(J, (k - 1L) * w + w)
    base[k, ix] <- 0.85
  }
  cl <- sample.int(K, n, TRUE)
  Y <- matrix(0L, n, J)
  for (j in seq_len(J)) Y[, j] <- 1L + stats::rbinom(n, 1L, base[cl, j])
  dd <- data.frame(x1 = stats::rnorm(n)); dd$Y <- Y
  list(dd = dd, Y = Y)
}
grid <- expand.grid(K = c(2L, 3L), n = c(200L, 600L),
                    design = c("mono", "block"), stringsAsFactors = FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  do.call(rbind, lapply(1:3, function(r) {
    K <- grid$K[i]; n <- grid$n[i]
    seed <- 31000L + 100L * i + r
    s <- if (grid$design[i] == "mono") sim_mono(seed, n, K) else
      sim_block(seed, n, K)
    nc <- rep(2L, ncol(s$Y))
    ll <- function(st) {
      f <- try(suppressWarnings(suppressMessages(
        if (is.null(st)) frm(bf(Y ~ x1), family = lca(K = K),
                             data = s$dd, control = tight)
        else frm(bf(Y ~ x1), family = lca(K = K), data = s$dd,
                 control = tight, start = st))), silent = TRUE)
      if (inherits(f, "try-error")) NA_real_ else as.numeric(logLik(f))
    }
    new <- ll(NULL)
    old <- ll(profiles_from(s$Y, nc, K, slice_score(s$Y, nc, K), w = 0))
    # a 10-refit perturbed-start multistart as the local reference
    set.seed(seed + 7L)
    f0 <- suppressWarnings(suppressMessages(
      frm(bf(Y ~ x1), family = lca(K = K), data = s$dd, control = tight)))
    p0 <- f0$frame$par_template[f0$frame$extra_names]
    ms <- max(vapply(1:10, function(z)
      ll(lapply(p0, function(v) v + stats::rnorm(length(v)))), numeric(1)),
      na.rm = TRUE)
    data.frame(design = grid$design[i], K = K, n = n, seed = seed,
               ll_new = new, ll_old = old, ll_best10 = ms)
  }))
}))
res$best <- pmax(res$ll_new, res$ll_old, res$ll_best10, na.rm = TRUE)
res$new_lost <- res$ll_new < res$best - 1e-6 * abs(res$best)
res$old_lost <- res$ll_old < res$best - 1e-6 * abs(res$best)
print(res, row.names = FALSE, digits = 9)
cat("\nreplicates:", nrow(res), "\n")
cat("the NEW start lost:", sum(res$new_lost), "\n")
cat("the OLD start lost:", sum(res$old_lost), "\n")
cat("cases where the new start is WORSE than the old:",
    sum(res$ll_new < res$ll_old - 1e-6 * abs(res$ll_old), na.rm = TRUE),
    "\n")
cat("cases where the new start is BETTER than the old:",
    sum(res$ll_new > res$ll_old + 1e-6 * abs(res$ll_old), na.rm = TRUE),
    "\n")
if (any(res$new_lost)) {
  cat("\nwhere the new start lost:\n")
  print(res[res$new_lost, ], row.names = FALSE, digits = 10)
}
