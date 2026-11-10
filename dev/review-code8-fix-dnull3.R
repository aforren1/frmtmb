# Fix-round check: headroom_ritz_flat() now also takes se_design_null()
# directions, along which the objective is asked only by straight probes.
# A design null is flat in the LIKELIHOOD, not in a penalized objective:
# y ~ x + s(x) puts s(x).fx1 on x's exact alias, and a prior on the two
# coefficients curves that direction. A stop short along it (data
# directions at their maximum) is made by an optimizer that runs nlminb
# to the MAP and then steps delta along the null direction, reported
# with code 8 or code 0.
#   RV_ARM=before|after Rscript dev/review-code8-fix-dnull.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
set.seed(5)
n <- 200
d <- data.frame(x = runif(n, 0, 10))
d$y <- 1 + 0.5 * d$x + sin(d$x) + rnorm(n, 0, 0.5)
form <- bf(y ~ x + s(x))
rows <- list()
for (psd in c(1, 10, 100, 1000, 10000)) {
  pr <- set_prior(sprintf("normal(0, %g)", psd), class = "b", coef = "x") +
    set_prior(sprintf("normal(0, %g)", psd), class = "b", coef = "s(x).fx1")
  ref <- rfit(form, data = d, prior = pr)
  f <- ref$f
  nm <- ns$outer_par_names(f)
  N <- ns$se_design_null(f)
  cat("psd", psd, "ref conds:", cond_txt(ref), "\n  design null cols:",
      NCOL(N), " loads:", paste(sprintf("%s %.3g", nm[abs(N[, 1]) > 1e-8],
                                         N[abs(N[, 1]) > 1e-8, 1]),
                               collapse = ", "), "\n")
  if (!NCOL(N)) next
  v <- N[, 1]
  for (code in c(1L, 0L)) for (delta in c(0.5, 2, 5, 20, 50) * psd / 10) {
    opt <- function(par, fn, gr, lower, upper, control) {
      o <- stats::nlminb(par, fn, gr, control = list(iter.max = 1000,
                                                     eval.max = 2000))
      p <- o$par + delta * v
      list(par = p, objective = fn(p), convergence = code,
           message = if (code == 1L) "false convergence (8)" else
             "relative convergence (4)")
    }
    r <- rfit(form, data = d, prior = pr,
              control = frmtmb_control(optimizer = opt, restarts = 0))
    if (r$err) { cat("err", conditionMessage(r$f), "\n"); next }
    g <- r$f
    txt <- c(r$w, r$m)
    rows[[length(rows) + 1L]] <- data.frame(psd = psd, code = code,
      delta = delta, short = signif(g$opt$objective - f$opt$objective, 3),
      says_max = any(grepl("so the estimates are at the maximum", txt)),
      grad_w = any(grepl("Large maximum", r$w)),
      code_w = any(r$wc == "frmtmb_conv_code"),
      at_max_m = any(r$mc == "frmtmb_converged_at_max"),
      nwarn = length(r$w), conds = gsub("[	
]", " ", paste(c(paste0("W[", r$wc, "] ", r$w), paste0("M[", r$mc, "] ", r$m)), collapse = " || ")))
  }
}
X <- do.call(rbind, rows)
options(width = 240)
print(X, right = FALSE)
write.table(X, sprintf("dev/review-code8-out/fix-dnull3-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("short > 1e-3:", sum(X$short > 1e-3), "; called at max:",
    sum(X$short > 1e-3 & (X$says_max | X$at_max_m)),
    "; no short verdict (no gradient warning, not code-warned as short):",
    sum(X$short > 1e-3 & !X$grad_w & (X$says_max | !X$code_w)), "\n")
cat("DONE\n")
