# Review of lane arspat2, attack 3: cost of the margin route at n = 1000
# (25 x 40 rook lattice, row-standardized), m = 50 and 200, (1 | h) 8
# levels, student() with gamma(2, 0.1) on nu. One fit per cell (seconds,
# elapsed), the tape's size (obj$env$...), and the R process's peak
# working set (ps). Base when REVLIB="". Arg: "fcor", "sar" or both.
source("dev/reviews/arspat2-review-common.R")
args <- commandArgs(TRUE)
which_t <- if (length(args)) args else c("sar", "fcor")
ms <- as.integer(strsplit(Sys.getenv("REV_MS", "50,200"), ",")[[1]])
g <- expand.grid(r = seq_len(25), c = seq_len(40))
n <- nrow(g)
Aj <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(n, n))
Aj <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- Aj / rowSums(Aj)
set.seed(11)
d <- data.frame(x = rnorm(n), h = factor(rep(1:8, length.out = n)))
A0 <- diag(n) - 0.4 * W
d$y <- as.vector(solve(A0, 1 + 0.5 * d$x + rnorm(8)[d$h] +
                         rnorm(n) * sqrt(4 / rchisq(1, 4))))
M <- if ("fcor" %in% which_t) solve(crossprod(A0)) else NULL
pr <- set_prior("gamma(2, 0.1)", class = "nu")
peak <- function() {
  if (requireNamespace("ps", quietly = TRUE)) {
    mi <- ps::ps_memory_info()
    v <- mi[["peak_wset"]] %||% mi[["rss"]]
    sprintf("%.2f GB", v / 2^30)
  } else "NA"
}
for (m in ms) {
  dm <- d
  set.seed(m)
  dm$y[sort(sample(n, m))] <- NA
  for (ty in which_t) {
    fo <- if (ty == "fcor") {
      y | mi() ~ x + (1 | h) + fcor(M)
    } else {
      y | mi() ~ x + (1 | h) + sar(W, type = "error")
    }
    d2 <- if (ty == "fcor") list(M = M) else list(W = W)
    gc()
    t0 <- proc.time()
    fit <- tryCatch(quiet(frm(fo, data = dm, data2 = d2, family = student(),
                              prior = pr)), error = function(e) e)
    tt <- proc.time() - t0
    if (inherits(fit, "error")) {
      cat(sprintf("n %d m %d %s: ERROR %s\n", n, m, ty, conditionMessage(fit)))
      next
    }
    # one objective + gradient evaluation at the optimum, timed alone
    t1 <- proc.time()[["elapsed"]]
    fit$obj$fn(fit$opt$par)
    fit$obj$gr(fit$opt$par)
    t_eval <- proc.time()[["elapsed"]] - t1
    cat(sprintf("n %d m %d %s: fit %.2f s elapsed (%.2f user), evals %s, fn+gr %.3f s, logLik %.4f, conv %d, peak %s\n",
                n, m, ty, tt[["elapsed"]], tt[["user.self"]],
                paste(fit$opt$evaluations, collapse = "/"), t_eval,
                as.numeric(logLik(fit)), fit$opt$convergence, peak()))
    rm(fit)
  }
}
