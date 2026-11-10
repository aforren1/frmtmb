# Review of lane code8, item 1: stops on a curved valley cut by optim()'s
# BFGS at maxit = k and reported as nlminb's "false convergence (8)", as
# the lane's own test does with nlminb (test-se-conv.R, conv_fake8()).
# Prints every condition in full for the cases that matter.
#   RV_ARM=before|after Rscript dev/review-code8-valley8.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
mk8 <- function(k, method = "BFGS") {
  function(par, fn, gr, lower, upper, control) {
    o <- stats::optim(par, fn, gr, method = method,
                      control = list(maxit = k))
    list(par = o$par, objective = o$value, convergence = 1L,
         message = "false convergence (8)")
  }
}
res <- list()
shown <- 0L
for (x0 in c(30, 100)) for (seed in 1:6) {
  set.seed(seed)
  n <- 60
  d <- data.frame(x = runif(n, x0, x0 + 1))
  d$y <- 2 * exp((3 / x0) * d$x) + rnorm(n, 0, 1)
  form <- bf(y ~ A * exp(B * x), A ~ 1, B ~ 1, nl = TRUE)
  ref <- rfit(form, data = d)
  fref <- ref$f$opt$objective
  for (k in c(seq(10, 100, 5), 150, 200, 300, 500, 1000)) {
    r <- rfit(form, data = d,
              control = frmtmb_control(optimizer = mk8(k), restarts = 0))
    if (r$err) next
    f <- r$f
    short <- f$opt$objective - fref
    txt <- c(r$w, r$m)
    room <- tryCatch(ns$conv_headroom(f, drop(f$obj$gr(f$opt$par))),
                     error = function(e) NA)
    s <- se_of(f)
    row <- data.frame(x0 = x0, seed = seed, k = k, short = signif(short, 4),
      room = signif(room, 3),
      says_max = any(grepl("so the estimates are at the maximum", txt)),
      grad_w = any(grepl("Large maximum absolute gradient", txt)),
      at_max_m = any(r$mc == "frmtmb_converged_at_max"),
      nwarn = length(r$w), se_lost = length(s$lost),
      se_B = signif(s$se[["B_(Intercept)"]], 4))
    res[[length(res) + 1L]] <- row
    if ((row$at_max_m || row$says_max) && short > 1e-3 && shown < 4) {
      shown <- shown + 1L
      cat("---- x0", x0, "seed", seed, "k", k, "short", short, "\n")
      cat(paste0("  W: ", r$w, collapse = "\n"), "\n")
      cat(paste0("  M[", r$mc, "]: ", r$m, collapse = "\n"), "\n")
    }
  }
}
X <- do.call(rbind, res)
options(width = 220)
print(X, right = FALSE)
write.table(X, sprintf("dev/review-code8-out/valley8-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("short > 2e-5 with the converged_at_max message:",
    sum(X$at_max_m & X$short > 2e-5), "of", sum(X$short > 2e-5), "\n")
cat("short > 1e-3 whose warning says 'at the maximum':",
    sum(X$says_max & X$short > 1e-3), "\n")
cat("short > 1e-3 with no warning:", sum(X$nwarn == 0 & X$short > 1e-3), "\n")
cat("DONE\n")
