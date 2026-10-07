# Reviewer of lane setier: the c0k curved ridge over seeds. Per seed,
# which of a_(Intercept), k_(Intercept), theta_1 (unidentified: the
# Laplace likelihood is invariant to (a0, k, b) -> (s a0, k / s, s b))
# keep a finite SE, and the warnings.
#   Rscript dev/setier-rev-c0k-sweep.R lane|base [seeds]
args <- commandArgs(TRUE)
arm <- args[1]
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:20
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, find.package("frmtmb"), "\n")
X <- list()
for (s in seeds) {
  set.seed(s)
  G <- 8
  dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
  dn$x <- rnorm(nrow(dn))
  u <- rnorm(G, 0, 0.6)[dn$g]
  dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(
    frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
           nl = TRUE), data = dn),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(x) invokeRestart("muffleMessage")),
    error = function(e) NULL)
  if (is.null(f)) { cat("seed", s, "error\n"); next }
  nm <- ns$outer_par_names(f)
  se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
  names(se) <- nm
  p <- f$opt$par
  f0 <- f$obj$fn(p)
  ia <- match("a_(Intercept)", nm); ik <- match("k_(Intercept)", nm)
  it <- match("theta_1", nm)
  q <- p; q[ia] <- 2 * p[ia]; q[ik] <- p[ik] / 2; q[it] <- p[it] + log(2)
  ridge <- f$obj$fn(q) - f0
  rid <- c("a_(Intercept)", "k_(Intercept)", "theta_1")
  cat(sprintf("seed %2d code %d | ridge change %.2g | SE a %.3g k %.3g theta_1 %.3g | %s\n",
              s, f$opt$convergence, ridge, se[[rid[1]]], se[[rid[2]]],
              se[[rid[3]]], if (length(w)) substr(w[1], 1, 60) else "silent"))
  X[[length(X) + 1L]] <- data.frame(seed = s, code = f$opt$convergence,
    flat = abs(ridge) < 1e-6, nfin = sum(is.finite(se[rid])),
    warned = length(w) > 0)
}
X <- do.call(rbind, X)
x <- X[X$flat & X$code == 0, ]
cat(sprintf("code-0 fits with a flat ridge: %d | with a finite SE on a ridge parameter: %d | of those silent: %d\n",
            nrow(x), sum(x$nfin > 0), sum(x$nfin > 0 & !x$warned)))
