source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# A degenerate fit (MLE at infinity) under an engaged default: does any
# warning survive? Seeds 511..520: a poisson with one all-zero factor
# cell, and a bernoulli with complete separation on a second covariate,
# each with a covariate x at sd 1e-4 so the default engages.
run <- function(...) {
  w <- character(0)
  f <- withCallingHandlers(tryCatch(frm(...), error = function(e) e),
         warning = function(x) { w <<- c(w, conditionMessage(x))
                                 invokeRestart("muffleWarning") })
  list(n = length(w), w = unique(substr(sub("[(].*", "", w), 1, 60)),
       code = if (inherits(f, "error")) NA else f$opt$convergence)
}
tab <- NULL
for (seed in 511:520) {
  set.seed(seed)
  n <- 240
  f <- factor(rep(c("a", "b", "c"), each = n / 3))
  x <- rnorm(n) * 1e-4
  yp <- rpois(n, c(a = 3, b = 2, c = 0)[as.character(f)] * exp(2000 * x))
  z <- rnorm(n); yb <- as.integer(z > 0)
  d <- data.frame(f, x, yp, yb, z)
  for (case in c("poisson zero cell", "bernoulli separated")) {
    fo <- if (case == "poisson zero cell") yp ~ f + x else yb ~ z + x
    fam <- if (case == "poisson zero cell") poisson() else bernoulli()
    a <- run(fo, family = fam, data = d)
    b <- run(fo, family = fam, data = d,
             control = frmtmb_control(autoscale = FALSE))
    tab <- rbind(tab, data.frame(seed, case, warnD = a$n, warnF = b$n,
                                 codeD = a$code, codeF = b$code,
                                 msgD = paste(a$w, collapse = "; "),
                                 msgF = paste(b$w, collapse = "; ")))
  }
}
options(width = 250)
print(tab, right = FALSE)
cat("\nfits with >= 1 warning under FALSE but 0 under default:",
    sum(tab$warnF > 0 & tab$warnD == 0), "of", nrow(tab), "\n")
