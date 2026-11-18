# Lane nbrobust: the cost of the negbinomial density.
#   Rscript dev/nbrobust-cost.R dens          one process, the two densities
#   Rscript dev/nbrobust-cost.R glmm <base|lib>   a negbinomial GLMM fit
# dens: RTMB::dnbinom_robust() against frmtmb:::nb_lpmf() on a taped
#   -sum(log mass) of 2000 rows, mean from a linear predictor, one shape
#   ("scalar") or a shape per row ("vector"), interleaved in blocks of at
#   least 1.2 s, the minimum over 5 rounds; the control is the robust arm
#   timed twice, and must read 1.
# glmm: y ~ x + (1 | g), 2000 rows, 50 groups, shape 2 (seed 1); the
#   fit's obj$fn and obj$gr at the optimum, and the whole frm() call,
#   with an arithmetic control in the same process.
args <- commandArgs(TRUE)
mode <- args[1]
arm <- if (length(args) >= 2) args[2] else "lib"
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-nbrobust-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
block <- function(f, min_t = 1.2) {
  reps <- 1L
  repeat {
    t <- system.time(for (i in seq_len(reps)) f())[["elapsed"]]
    if (t >= min_t) return(t / reps)
    reps <- reps * 2L
  }
}
control <- function() {
  x <- 0
  for (i in 1:2e5) x <- x + sqrt(i)
  x
}
set.seed(1)
n <- 2000
g <- factor(sample(sprintf("g%02d", 1:50), n, TRUE))
x <- stats::rnorm(n)
u <- stats::rnorm(50, 0, 0.5)
mu <- exp(1 + 0.4 * x + u[as.integer(g)])
y <- stats::rnbinom(n, size = 2, mu = mu)
cat("counts: max", max(y), " mean", mean(y), " rows above 24", sum(y > 24),
    "\n")
if (mode == "dens") {
  X <- cbind(1, x)
  mk <- function(dens, vec) {
    f <- function(p) {
      lmu <- as.vector(X %*% p$b)
      lsh <- if (vec) p$ls + 0 * lmu else p$ls
      -sum(dens(y, lmu, lsh))
    }
    RTMB::MakeADFun(f, list(b = c(1, 0.4), ls = log(2)), silent = TRUE)
  }
  robust <- function(y, lmu, lsh) {
    RTMB::dnbinom_robust(y, lmu, 2 * lmu - lsh, log = TRUE)
  }
  objs <- list(robust = mk(robust, FALSE), robust2 = mk(robust, FALSE),
               nb = mk(frmtmb:::nb_lpmf, FALSE),
               robust_v = mk(robust, TRUE), nb_v = mk(frmtmb:::nb_lpmf, TRUE))
  p <- objs$robust$par
  cat(sprintf("value check: robust %.12f nb %.12f\n", objs$robust$fn(p),
              objs$nb$fn(p)))
  res <- list()
  for (round in 1:5) {
    for (nm in names(objs)) {
      o <- objs[[nm]]
      res[[length(res) + 1L]] <- data.frame(
        round = round, arm = nm, fn = block(function() o$fn(p)),
        gr = block(function() o$gr(p)))
    }
  }
  r <- do.call(rbind, res)
  m <- aggregate(cbind(fn, gr) ~ arm, r, min)
  rownames(m) <- m$arm
  for (nm in rownames(m)) {
    den <- if (grepl("_v$", nm)) "robust_v" else "robust"
    cat(sprintf("%-9s fn %.3e s gr %.3e s  ratio to %s: fn %.3f gr %.3f\n",
                nm, m[nm, "fn"], m[nm, "gr"], den,
                m[nm, "fn"] / m[den, "fn"], m[nm, "gr"] / m[den, "gr"]))
  }
}
if (mode == "glmm") {
  d <- data.frame(y, x, g)
  tc <- block(control)
  tfit <- system.time(fit <- frm(y ~ x + (1 | g), family = negbinomial(),
                                 data = d))[["elapsed"]]
  p <- fit$opt$par
  cat(sprintf("logLik %.10f shape %.6f code %d\n", as.numeric(logLik(fit)),
              p[length(p)], fit$opt$convergence))
  res <- list()
  for (round in 1:5) {
    res[[round]] <- c(control = block(control),
                      fn = block(function() fit$obj$fn(p)),
                      gr = block(function() fit$obj$gr(p)))
  }
  r <- do.call(rbind, res)
  mins <- apply(r, 2, min)
  cat(sprintf(paste("RESULT glmm %s control %.4e fn %.4e gr %.4e",
                    "fn/control %.4f gr/control %.4f fit %.2f\n"), arm,
              mins[["control"]], mins[["fn"]], mins[["gr"]],
              mins[["fn"]] / mins[["control"]],
              mins[["gr"]] / mins[["control"]], tfit))
}
cat("END\n")
