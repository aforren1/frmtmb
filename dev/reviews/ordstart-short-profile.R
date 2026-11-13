# Short check of lane ordstart's Q4 change: the five sweep fits whose
# unconverged step-function default the Q4 rule replaced by a converged
# proper start. At each new point, profile the log-likelihood along its
# largest coordinate (that coordinate fixed at t times its value, the
# rest re-optimized, warm-started outward from t = 1), and print the
# degeneracy table, the Hessian's smallest eigenvalues and the SE of that
# coordinate. An interior maximum falls off on both sides; a near-step
# function stays flat or rises as t grows.
#   Rscript ordstart-short-profile.R
.libPaths(c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/reviews/ordstart-review-sweep.R")
eval(parse(text = src[grep("^gen <- ", src):(grep("^cases <- list", src) - 1L)]))
q <- function(e) suppressWarnings(suppressMessages(e))
cases <- list(
  cs_sr_acat_seed4 = list(seed = 4, f = bf(y ~ cs(x) + z),
                          fam = mixture(sratio(), acat())),
  cs_sr_acat_seed1 = list(seed = 1, f = bf(y ~ cs(x) + z),
                          fam = mixture(sratio(), acat())),
  cs_sr_acat_seed3 = list(seed = 3, f = bf(y ~ cs(x) + z),
                          fam = mixture(sratio(), acat())),
  dpar_formula_seed4 = list(seed = 4, f = bf(y ~ 1, mu1 ~ x, mu2 ~ x + z),
                            fam = mixture(cumulative(), sratio())),
  dpar_formula_seed7 = list(seed = 7, f = bf(y ~ 1, mu1 ~ x, mu2 ~ x + z),
                            fam = mixture(cumulative(), sratio())))
ts <- c(0.25, 0.5, 0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.25, 1.5, 2, 3, 5, 10, 30)
op <- options(width = 170)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  d <- gen(cs$seed, n = 300, K = 4)
  f <- q(frm(cs$f, family = cs$fam, data = d))
  obj <- f$obj
  p <- f$opt$par
  ll0 <- as.numeric(logLik(f))
  vn <- tryCatch(variables(f), error = function(e) NULL)
  j <- which.max(abs(p))
  cat(sprintf("\n=== %s ll %.4f code %d max|par| %.4g at par %d (%s) %s\n",
              nm, ll0, f$opt$convergence, abs(p[j]), j, names(p)[j],
              if (length(vn) == length(p)) vn[j] else ""))
  print(frmtmb:::mixture_ord_degeneracy(f, "y"), digits = 4)
  H <- tryCatch(obj$he(p), error = function(e) NULL)
  if (!is.null(H)) {
    ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
    V <- tryCatch(solve((H + t(H)) / 2), error = function(e) NULL)
    cat("Hessian eigenvalues, smallest 3:", signif(tail(ev, 3), 3),
        "| largest:", signif(ev[1], 3), "| SE of that coordinate:",
        if (is.null(V)) NA else signif(sqrt(V[j, j]), 3), "\n")
  }
  lo <- f$cache$bounds$lower %||% rep(-Inf, length(p))
  hi <- f$cache$bounds$upper %||% rep(Inf, length(p))
  lo <- rep_len(lo, length(p))
  hi <- rep_len(hi, length(p))
  prof_at <- function(v, from) {
    fr <- function(qq) {
      pp <- p
      pp[-j] <- qq
      pp[j] <- v
      r <- obj$fn(pp)
      if (is.finite(r)) r else 1e10
    }
    gr <- function(qq) {
      pp <- p
      pp[-j] <- qq
      pp[j] <- v
      obj$gr(pp)[-j]
    }
    o <- q(nlminb(from, fr, gr, lower = lo[-j], upper = hi[-j],
                  control = list(iter.max = 2000, eval.max = 4000)))
    list(ll = -o$objective, par = o$par, code = o$convergence)
  }
  res <- numeric(length(ts))
  for (dir in c("up", "down")) {
    idx <- if (dir == "up") which(ts >= 1) else rev(which(ts <= 1))
    from <- p[-j]
    for (i in idx) {
      r <- prof_at(ts[i] * p[j], from)
      from <- r$par
      res[i] <- r$ll + f$opt$objective
    }
  }
  cat("profile logLik - logLik at t * value:\n")
  print(setNames(signif(res, 4), paste0("t", ts)))
}
cat("DONE\n")
