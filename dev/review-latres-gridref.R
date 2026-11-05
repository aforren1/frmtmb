# Is the exact reference exact? The forward likelihood at a coarse and a fine
# grid, at the true point and at the point the optimizer picked, per design.
root <- "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/"
for (design in c("gamma_T", "nb_T", "gamma_F", "probit_F", "pois_F", "bern_F")) {
  commandArgs <- function(trailingOnly = FALSE) c(design, "1", "1", "x")
  src <- readLines(paste0(root, "review-latres-sim2.R"))
  env <- new.env(parent = globalenv())
  eval(parse(text = src[seq_len(grep("^cols <- ", src) - 1L)]), env)
  D <- env$D
  for (seed in 1:3) {
    d <- env$simd(seed)
    st <- c(D$b0, D$b1, log(D$s), if (D$kind == "T") D$ar / sqrt(1 - D$ar^2) else qlogis((D$ar + 1) / 2), if (!is.null(D$shape)) log(D$shape))
    a <- env$exact_nll(st, d, ng = env$ngrid)
    b <- env$exact_nll(st, d, ng = 2 * env$ngrid + 1)
    cat(sprintf("%-9s seed %d  nll ng=%d %.6f  ng=%d %.6f  diff %.2e\n", design, seed, env$ngrid, a, 2 * env$ngrid + 1, b, a - b))
  }
}
commandArgs <- function(trailingOnly = FALSE) c("bern", "1", "1", "x")
src <- readLines(paste0(root, "review-latres-sim.R"))
env <- new.env(parent = globalenv())
eval(parse(text = src[seq_len(grep("^fit_one <- ", src) - 1L)]), env)
for (seed in 1:2) {
  d <- env$simd(seed)
  for (par in list(c(.2, .5, 0, .75), c(5, 1.7, log(10), .39), c(.27, .22, log(.15), 3.7))) {
    cat(sprintf("bern seed %d s=%.2f  nll 121 %.6f  301 %.6f\n", seed, exp(par[3]), env$exact_nll(par, d), env$exact_nll(par, d, ngrid = 301)))
  }
}
