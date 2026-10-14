# Reviewer of lane illse: degree-6 seeds 4 and 7 of dev/illse-rev-edge4.R
# (sd at zero): the warnings in full, the log-likelihood against lm(),
# and the theta SE.
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
for (seed in c(4, 7)) {
  deg <- 6
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 2), g = factor(rep(1:20, length.out = n)))
  d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "), "+ (1 | g)"))
  r <- cap(frm(fo, data = d))
  f <- r$v
  lmf <- stats::lm(stats::as.formula(paste("y ~", paste(xn, collapse = " + "))), data = d)
  cat(sprintf("seed %d | code %d | logLik(lm) - logLik(fit) %.2e | polish %s | explained %s\n", seed,
              f$opt$convergence, as.numeric(logLik(lmf)) - as.numeric(logLik(f)),
              paste(f$opt$polish, collapse = ","), format(f$cache$se_explained %||% "-")))
  cat(paste(" W:", r$w), sep = "\n")
  print(sqrt(diag(suppressWarnings(vcov(f, full = TRUE)))))
  print(f$opt$par)
}
