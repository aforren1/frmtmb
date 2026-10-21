# Reviewer: the runaway-phi claim on fresh seeds, both arms in one
# process (lgamma_shift_diff() and the plain lgamma difference brms
# writes, swapped with assignInNamespace before each fit is taped).
#   Rscript dev/nb2dm-rev-limit.R <first seed> <last seed> <K> <Nhi>
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
a <- commandArgs(trailingOnly = TRUE)
seeds <- seq(as.integer(a[1]), as.integer(a[2]))
K <- as.integer(a[3])
Nhi <- as.integer(a[4])
real <- frmtmb:::lgamma_shift_diff
plain <- function(a, s) lgamma(a + s) - lgamma(a)
fit_arm <- function(fun, d) {
  utils::assignInNamespace("lgamma_shift_diff", fun, "frmtmb")
  on.exit(utils::assignInNamespace("lgamma_shift_diff", real, "frmtmb"))
  msgs <- character()
  m <- withCallingHandlers(
    frm(Y | trials(N) ~ x, family = dirichlet_multinomial(), data = d),
    warning = function(w) {
      msgs <<- c(msgs, paste("W:", substr(conditionMessage(w), 1, 50)))
      invokeRestart("muffleWarning")
    },
    message = function(w) {
      msgs <<- c(msgs, paste("M:", substr(conditionMessage(w), 1, 50)))
      invokeRestart("muffleMessage")
    })
  list(m = m, msgs = paste(unique(msgs), collapse = " | "))
}
rows <- list()
for (seed in seeds) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = rnorm(n))
  E <- cbind(0, sapply(2:K, function(k) rnorm(1, 0, 0.5) +
                         rnorm(1, 0, 0.5) * d$x))
  P <- exp(E) / rowSums(exp(E))
  d$N <- sample(5:Nhi, n, TRUE)
  Y <- t(vapply(seq_len(n), function(i) {
    as.numeric(rmultinom(1, d$N[i], P[i, ]))
  }, numeric(K)))
  colnames(Y) <- letters[seq_len(K)]
  d$Y <- Y
  mn <- frm(Y | trials(N) ~ x, family = multinomial(), data = d)
  llm <- as.numeric(logLik(mn))
  shift_obj <- NULL
  for (arm in c("shift", "plain")) {
    r <- fit_arm(if (arm == "shift") real else plain, d)
    if (arm == "shift") shift_obj <- r$m$obj
    # each arm's estimate scored on the ACCURATE tape (the shift arm's):
    # what the plain arm's point is really worth
    ll_true <- -shift_obj$fn(r$m$opt$par)
    rows[[length(rows) + 1L]] <- data.frame(
      seed = seed, arm = arm, excess = as.numeric(logLik(r$m)) - llm,
      excess_true = ll_true - llm,
      log_phi = r$m$estimates$betad[[1]], conv = r$m$opt$convergence,
      msg = r$m$opt$message, cond = r$msgs)
  }
}
out <- do.call(rbind, rows)
print(out[, 1:7], digits = 6, row.names = FALSE)
cat("conditions:\n")
print(unique(out$cond))
for (arm in c("shift", "plain")) {
  o <- out[out$arm == arm, ]
  run <- o$log_phi > 12
  cat(sprintf(paste("ARM %s K=%d Nhi=%d seeds %d-%d: run-off %d of %d;",
                    "of those excess > 1e-8: %d, max excess %.3g,",
                    "conv != 0: %d; finite-phi conv != 0: %d\n"),
              arm, K, Nhi, min(seeds), max(seeds), sum(run), nrow(o),
              sum(o$excess[run] > 1e-8),
              if (any(run)) max(o$excess[run]) else NA,
              sum(o$conv[run] != 0), sum(o$conv[!run] != 0)))
}
