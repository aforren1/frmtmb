# Fix-round check (m1): se_check_at_fit() budgets 4 np gradients. Which
# fits with more than 10 outer parameters now defer the SE check, and does
# a verdict move from frm() to summary()? Models with a known SE verdict
# (a slope sd at 0: the boundary message; an exact nonlinear alias: the
# SE warning), 11 to 14 outer parameters, counted evaluations printed.
#   RV_ARM=before|after Rscript dev/review-code8-fix-defer2.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
catch <- function(expr) {
  w <- character(0); m <- character(0)
  v <- withCallingHandlers(expr,
    warning = function(x) { w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning") },
    message = function(x) { m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage") })
  list(v = v, txt = substr(gsub("\n", " ", c(paste("W:", w), paste("M:", m))),
                           1, 70))
}
for (seed in 5:60) for (kx in c(6)) {
  set.seed(seed)
  n <- 600
  d <- data.frame(g = factor(sample(1:30, n, TRUE)))
  for (j in 1:kx) d[[paste0("x", j)]] <- rnorm(n)
  d$y <- 1 + 0.3 * d$x1 + rnorm(30, 0, 0.5)[d$g] + rnorm(n)
  fo <- as.formula(paste("y ~", paste0("x", 1:kx, collapse = " + "),
                         "+ (1 + x1 | g)"))
  a <- catch(frm(fo, data = d))
  f <- a$v
  b <- catch(summary(f))
  cat(sprintf("seed %d kx %d np %d evals %s at_fit %s\n  frm: %s\n  summary: %s\n",
              seed, kx, length(f$opt$par), format(f$opt$evals),
              ns$se_check_at_fit(f),
              paste(a$txt, collapse = " | "), paste(b$txt, collapse = " | ")))
}
cat("DONE\n")
