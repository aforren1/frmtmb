# Lane selimit2, backlog 2754 and 2759: the c0k curved ridge
# (dev/seflags-rev2-c0k-sweep.R) rerun on 0.82.0 and the lane. The
# Laplace likelihood is invariant to (a0, k, b) -> (s a0, k / s, s b),
# so a_(Intercept), k_(Intercept) and theta_1 have no standard error.
# 2754: seeds 1 and 7 (17 under OpenBLAS) kept one of the three finite.
# 2759: seeds 3 and 7 kept zeroed loadings from an unconfirmed downward
# direction. Per seed: the ridge's change in the objective, the three
# SEs, the lost reasons and the warnings.
#   Rscript dev/selimit2-c0k.R <base|lane> [seeds]
args <- commandArgs(TRUE)
arm <- args[1]
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:20
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, dirname(find.package("frmtmb")), "| BLAS",
    basename(extSoftVersion()[["BLAS"]]), La_library(), "\n")
X <- list()
for (s in seeds) {
  set.seed(s)
  G <- 8
  dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
  dn$x <- rnorm(nrow(dn))
  u <- rnorm(G, 0, 0.6)[dn$g]
  dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
  w <- character()
  m <- character()
  f <- tryCatch(withCallingHandlers(
    frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
           nl = TRUE), data = dn),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(x) {
      m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
    }),
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
  lost <- ns$sdr_of(f)$se_lost
  cat(sprintf(paste0("seed %2d code %d | ridge change %.2g | SE a %.3g k %.3g ",
                     "theta_1 %.3g | lost %s | warn %d msg %d | %s\n"),
              s, f$opt$convergence, ridge, se[[rid[1]]], se[[rid[2]]],
              se[[rid[3]]],
              if (length(lost)) paste(names(lost), lost, sep = ":",
                                      collapse = ",") else "none",
              length(w), length(m),
              if (length(w)) substr(w[1], 1, 70) else "silent"))
  X[[length(X) + 1L]] <- data.frame(seed = s, code = f$opt$convergence,
    flat = abs(ridge) < 1e-6, nfin = sum(is.finite(se[rid])),
    warned = length(w) > 0)
}
X <- do.call(rbind, X)
x <- X[X$flat, ]
cat(sprintf(paste0("fits with a flat ridge: %d | with a finite SE on a ridge ",
                   "parameter: %d (seeds %s) | of those silent: %d\n"),
            nrow(x), sum(x$nfin > 0),
            paste(x$seed[x$nfin > 0], collapse = ","),
            sum(x$nfin > 0 & !x$warned)))
