# Review of lane code8, items 1 and 2, with real stops: y = A exp(B x),
# x on [30, 31], fits cut short by nlminb's own controls (iter.max,
# rel.tol), each compared with the default fit. What does each build say?
#   RV_ARM=before|after Rscript dev/review-code8-valley-fit.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
res <- list()
for (x0 in c(30, 100)) for (seed in 1:3) {
  set.seed(seed)
  n <- 60
  d <- data.frame(x = runif(n, x0, x0 + 1))
  d$y <- 2 * exp((3 / x0) * d$x) + rnorm(n, 0, 1)
  form <- bf(y ~ A * exp(B * x), A ~ 1, B ~ 1, nl = TRUE)
  ref <- rfit(form, data = d)
  fref <- ref$f$opt$objective
  ctl <- c(lapply(c(2:12, 15, 20, 30), function(k) list(iter.max = k)),
           lapply(c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6), function(t) list(rel.tol = t)),
           lapply(c(1e-2, 1e-3, 1e-4), function(t) list(x.tol = t)))
  for (cc in ctl) {
    r <- rfit(form, data = d,
              control = frmtmb_control(optCtrl = c(cc, eval.max = 1000),
                                       restarts = 0))
    if (r$err) next
    f <- r$f
    short <- f$opt$objective - fref
    txt <- c(r$w, r$m)
    res[[length(res) + 1L]] <- data.frame(x0 = x0, seed = seed,
      ctl = paste(names(cc), cc[[1]]), code = f$opt$convergence,
      msg = substr(f$opt$message, 1, 22), short = signif(short, 4),
      says_max = any(grepl("so the estimates are at the maximum", txt)),
      says_short = any(grepl("short of the|would still gain|unusable", txt)),
      grad_w = any(grepl("Large maximum absolute gradient", txt)),
      code_w = any(r$wc == "frmtmb_conv_code"),
      at_max_m = any(r$mc == "frmtmb_converged_at_max"),
      nwarn = length(r$w), se_lost = length(se_of(f)$lost))
  }
}
X <- do.call(rbind, res)
options(width = 220)
print(X, right = FALSE)
write.table(X, sprintf("dev/review-code8-out/valley-fit-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
bad <- X$short > 1e-3 & !X$says_short & !X$grad_w
cat("short by > 1e-3 with no short/gradient verdict:", sum(bad), "of",
    sum(X$short > 1e-3), "\n")
cat("  of those silent (no warning at all):", sum(bad & X$nwarn == 0), "\n")
cat("  of those calling the stop the maximum:", sum(bad & X$says_max), "\n")
cat("  shortfalls:", paste(signif(X$short[bad], 3), collapse = " "), "\n")
cat("DONE\n")
