# dirichlet_multinomial() on data with no overdispersion, where phi runs
# to infinity and the model's supremum is the multinomial fit's
# log-likelihood. A fit that reports more than that has fitted noise in
# the lgamma differences (the student-t story of lgamma_shift_diff()).
# Also the same on beta_binomial(), the K = 2 case, as a control.
#   Rscript dev/nb2dm-dm-limit.R [lib] > dev/nb2dm-out/dm-limit.txt 2>&1
a <- commandArgs(trailingOnly = TRUE)
lib <- if (length(a)) a[1] else "C:/Users/adf44/source/r/wt-nb2dm-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r10",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
res <- list()
for (seed in 1:12) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = rnorm(n))
  E <- cbind(0, 0.3 + 0.5 * d$x, -0.2 - 0.4 * d$x)
  P <- exp(E) / rowSums(exp(E))
  d$N <- sample(5:40, n, TRUE)
  Y <- t(vapply(seq_len(n), function(i) {
    as.numeric(rmultinom(1, d$N[i], P[i, ]))
  }, numeric(3)))
  colnames(Y) <- c("a", "b", "c")
  d$Y <- Y
  mn <- frm(Y | trials(N) ~ x, family = multinomial(), data = d)
  w <- character()
  dm <- withCallingHandlers(
    frm(Y | trials(N) ~ x, family = dirichlet_multinomial(), data = d),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    },
    message = function(cnd) invokeRestart("muffleMessage"))
  lphi <- dm$estimates$betad[[1]]
  res[[seed]] <- data.frame(
    seed = seed, ll_mn = as.numeric(logLik(mn)),
    ll_dm = as.numeric(logLik(dm)),
    excess = as.numeric(logLik(dm)) - as.numeric(logLik(mn)),
    log_phi = lphi, conv = dm$opt$convergence,
    max_coef_gap = max(abs(unlist(fixef_by_dpar(dm)[c("mub", "muc")]) -
                             unlist(fixef_by_dpar(mn)))),
    warn = paste(substr(w, 1, 60), collapse = " | "))
}
out <- do.call(rbind, res)
print(out[, 1:7], digits = 10, row.names = FALSE)
cat("warnings:\n")
print(out$warn)
# a fit whose phi stays finite found real (chance) overdispersion and
# may beat the multinomial; one whose phi ran off is at the multinomial
# limit, which it can approach but not exceed
run <- out$log_phi > 12
cat("fits with phi run off (log phi > 12):", sum(run), "of", nrow(out),
    "\n")
cat("  of those, above the multinomial by more than 1e-8:",
    sum(out$excess[run] > 1e-8), "; max excess",
    format(max(out$excess[run]), digits = 6), "\n")
cat("  of those, optimizer code not 0:", sum(out$conv[run] != 0), "\n")
cat("fits with finite phi: optimizer code not 0:",
    sum(out$conv[!run] != 0), "\n")
