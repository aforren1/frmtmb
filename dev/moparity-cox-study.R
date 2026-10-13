# cox()'s baseline simplex at the maximum likelihood. The construction
# of dev/optima-rev-cox1.R (n = 300, right censoring), seeds 1..N. Per
# seed: the fit's log-likelihood, its smallest weight, the Frank-Wolfe
# gap of the baseline simplex at the fit's coefficients (the most the
# log-likelihood can still rise by moving the weights alone, an upper
# bound because it is concave in them), and the log-likelihood of a
# polished maximum: alternating exponentiated-gradient steps on the
# weights (to a gap below 1e-11) and Newton steps on the coefficients,
# from the fit, until neither moves.
#   Rscript dev/moparity-cox-study.R lane|base N
args <- commandArgs(TRUE)
arm <- args[1]
N <- as.integer(args[2] %||% 200)
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
loglik_hand <- function(beta, s, X, Zb, Zc, ev) {
  eta <- as.vector(X %*% beta)
  sum((log(as.vector(Zb %*% s)) + eta)[ev]) -
    sum(as.vector(Zc %*% s) * exp(eta))
}
polish <- function(beta, s, X, Zb, Zc, ev) {
  # The intercept and the scale of s are one parameter (h0 exp(b0)), so
  # the intercept goes into s: s >= 0 with no sum constraint, maximized
  # by the multiplicative (EM) update for sum log(a_i s) - c s, which
  # rises at every step; the slopes by Newton, concave in eta.
  s <- pmax(s, 1e-300) * exp(beta[1])
  b <- beta[-1]
  Xs <- X[, -1, drop = FALSE]
  Zbe <- Zb[ev, , drop = FALSE]
  num_ev <- colSums(Xs[ev, , drop = FALSE])
  for (outer in 1:500) {
    mu <- exp(as.vector(Xs %*% b))
    cvec <- colSums(Zc * mu)
    for (it in 1:2000) {
      s <- s * colSums(Zbe / as.vector(Zbe %*% s)) / cvec
    }
    b_old <- b
    for (it in 1:50) {
      cb <- as.vector(Zc %*% s)
      mu <- exp(as.vector(Xs %*% b))
      gr <- num_ev - colSums(Xs * (cb * mu))
      H <- -crossprod(Xs * sqrt(cb * mu))
      step <- solve(H, gr)
      b <- b - step
      if (max(abs(step)) < 1e-14) break
    }
    if (max(abs(b - b_old)) < 1e-13) break
  }
  beta <- c(log(sum(s)), b)
  s <- s / sum(s)
  mu <- exp(as.vector(X %*% beta))
  g <- colSums(Zbe / as.vector(Zbe %*% s)) - colSums(Zc * mu)
  list(beta = beta, s = s, ll = loglik_hand(beta, s, X, Zb, Zc, ev),
       gap = max(g) - sum(s * g))
}
rows <- list()
for (seed in seq_len(N)) {
  set.seed(seed)
  n <- 300
  x <- rnorm(n)
  tt <- stats::rexp(n, exp(-0.5 + 0.7 * x))
  ct <- stats::rexp(n, 0.2)
  dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
  w <- character()
  fit <- withCallingHandlers(
    frm(bf(time | cens(cens) ~ x), family = cox(), data = dd),
    warning = function(z) {
      w <<- c(w, conditionMessage(z))
      invokeRestart("muffleWarning")
    })
  av <- fit$frame$aterm_values[[1L]]
  Zb <- av$Zbhaz
  Zc <- av$Zcbhaz
  X <- as.matrix(fit$frame$linpreds[[1L]]$X)
  ev <- dd$cens == 0
  s <- as.numeric(cox_baseline(fit))
  beta <- as.numeric(fit$estimates$beta)
  ll_fit <- as.numeric(logLik(fit))
  ll_hand <- loglik_hand(beta, s, X, Zb, Zc, ev)
  eta <- as.vector(X %*% beta)
  g <- colSums(Zb[ev, , drop = FALSE] / as.vector(Zb[ev, , drop = FALSE] %*%
                                                    s)) -
    colSums(Zc * exp(eta))
  fw <- max(g) - sum(s * g)
  p <- polish(beta, s, X, Zb, Zc, ev)
  p2 <- polish(beta, rep(1 / length(s), length(s)), X, Zb, Zc, ev)
  best <- max(p$ll, p2$ll)
  rows[[seed]] <- data.frame(
    seed = seed, code = fit$opt$convergence, ll = ll_fit,
    hand_gap = ll_fit - ll_hand, minw = min(s), fw = fw,
    gap = best - ll_fit, polish_minw = min(if (p$ll >= p2$ll) p$s else p2$s),
    polish_fw = if (p$ll >= p2$ll) p$gap else p2$gap,
    evals = sum(unlist(fit$opt$evals %||% NA)),
    warn = substr(paste(w, collapse = " || "), 1, 80))
  cat(sprintf("seed %3d code %d ll %.8f minw %.3g fw %.3g gap %.3g pfw %.2g warn %s\n",
              seed, fit$opt$convergence, ll_fit, min(s), fw, best - ll_fit,
              rows[[seed]]$polish_fw,
              rows[[seed]]$warn))
}
res <- do.call(rbind, rows)
saveRDS(res, file.path("dev/moparity-log", paste0("cox-study-", arm, ".rds")))
cat("\nSUMMARY", arm, ": seeds", nrow(res),
    "| gap > 1e-6:", sum(res$gap > 1e-6), "| > 1e-3:", sum(res$gap > 1e-3),
    "| max gap", format(max(res$gap), digits = 4),
    "| min weight < 1e-6 at the fit:", sum(res$minw < 1e-6),
    "| at the polished max:", sum(res$polish_minw < 1e-6),
    "| codes", paste(names(table(res$code)), table(res$code), sep = "=",
                     collapse = " "),
    "| any warning", sum(nzchar(res$warn)),
    "| max |logLik - hand|", format(max(abs(res$hand_gap)), digits = 3),
    "\n")
