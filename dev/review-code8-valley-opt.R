# Review of lane code8, item 2: natural short stops on a curved valley
# from optimizers a user can pass (optim()'s L-BFGS-B, BFGS, CG), each
# compared with the default nlminb fit. Which builds stay silent?
#   RV_ARM=before|after Rscript dev/review-code8-valley-opt.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
mkopt <- function(method) {
  function(par, fn, gr, lower, upper, control) {
    lo <- if (method == "L-BFGS-B") lower else -Inf
    hi <- if (method == "L-BFGS-B") upper else Inf
    o <- stats::optim(par, fn, gr, method = method, lower = lo, upper = hi,
                      control = list(maxit = 5000))
    list(par = o$par, objective = o$value,
         convergence = as.integer(o$convergence),
         message = o$message %||% "optim")
  }
}
res <- list()
for (x0 in c(10, 30, 100)) for (seed in 1:6) {
  set.seed(seed)
  n <- 60
  d <- data.frame(x = runif(n, x0, x0 + 1))
  d$y <- 2 * exp((3 / x0) * d$x) + rnorm(n, 0, 1)
  form <- bf(y ~ A * exp(B * x), A ~ 1, B ~ 1, nl = TRUE)
  ref <- rfit(form, data = d)
  fref <- ref$f$opt$objective
  for (meth in c("L-BFGS-B", "BFGS", "CG")) {
    r <- rfit(form, data = d,
              control = frmtmb_control(optimizer = mkopt(meth), restarts = 0))
    if (r$err) { cat(x0, seed, meth, "error", conditionMessage(r$f), "\n"); next }
    f <- r$f
    short <- f$opt$objective - fref
    txt <- c(r$w, r$m)
    res[[length(res) + 1L]] <- data.frame(x0 = x0, seed = seed, meth = meth,
      code = f$opt$convergence, short = signif(short, 4),
      says_max = any(grepl("so the estimates are at the maximum", txt)),
      grad_w = any(grepl("Large maximum absolute gradient", txt)),
      nwarn = length(r$w), se_lost = length(se_of(f)$lost),
      conds = substr(cond_txt(r), 1, 90))
  }
}
X <- do.call(rbind, res)
options(width = 220)
print(X, right = FALSE)
write.table(X, sprintf("dev/review-code8-out/valley-opt-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
bad <- X$short > 1e-3 & X$nwarn == 0
cat("short by > 1e-3 with no warning:", sum(bad), "of", sum(X$short > 1e-3),
    "\n  shortfalls:", paste(signif(X$short[bad], 3), collapse = " "), "\n")
cat("DONE\n")
