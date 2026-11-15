# Review of lane famend, item 3060: side effects of frechet_fit_check().
# Lane build only; the arm "off" replaces frechet_fit_check() by a no-op in
# this process, so the two arms differ only by the check.
#   Rscript dev/reviews/famend-review-frechet.R
lanelib <- "C:/Users/adf44/source/r/wt-famend-lib"
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(lanelib, base, ulib))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
real_check <- get("frechet_fit_check", ns)

cap <- function(expr) {
  w <- character(); m <- character()
  v <- withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
    })
  list(v = v, w = w, m = m)
}
dat <- function(seed, n = 300, nu0 = 1.1, G = 0) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / nu0) * rexp(n)^(-1 / nu0)
  if (G) {
    d$g <- factor(rep(seq_len(G), length.out = n))
    d$y <- d$y * exp(rnorm(G, 0, 0.3))[d$g]
  }
  d
}
state <- function(f) {
  env <- f$obj$env
  lf <- env$lfixed()
  list(best = env$last.par.best, last = env$last.par,
       vbest = env$value.best,
       bestfix = unname(env$last.par.best[lf]))
}
run <- function(label, fo, d, on, se = FALSE) {
  assignInNamespace("frechet_fit_check",
                    if (on) real_check else function(fit, resp) invisible(NULL),
                    ns)
  t0 <- proc.time()[["elapsed"]]
  r <- cap(frm(fo, family = frechet(), data = d, se = se))
  el <- proc.time()[["elapsed"]] - t0
  f <- r$v
  s <- state(f)
  V <- tryCatch(suppressWarnings(vcov(f, full = TRUE)), error = function(e) NULL)
  re <- tryCatch(unlist(ranef(f)), error = function(e) NULL)
  cat(sprintf("%-34s check %-3s elapsed %.2f s code %d objective %.12f\n",
              label, if (on) "on" else "off", el, f$opt$convergence,
              f$opt$objective))
  cat("   last.par.best (fixed) - opt$par:",
      paste(signif(s$bestfix - unname(f$opt$par), 3), collapse = " "), "\n")
  cat("   value.best - objective:", signif(s$vbest - f$opt$objective, 3), "\n")
  for (x in r$w) cat("   W:", substr(gsub("\\s+", " ", x), 1, 150), "\n")
  for (x in r$m) cat("   M:", substr(gsub("\\s+", " ", x), 1, 150), "\n")
  list(f = f, s = s, V = V, re = re, ll = as.numeric(logLik(f)),
       fit_pred = tryCatch(fitted(f)[1:3, 1], error = function(e) NA))
}
cmp <- function(a, b) {
  same <- function(x, y) identical(x, y)
  cat("   identical last.par.best:", same(a$s$best, b$s$best),
      " last.par:", same(a$s$last, b$s$last),
      " vcov:", same(a$V, b$V),
      " ranef:", same(a$re, b$re),
      " logLik:", same(a$ll, b$ll),
      " fitted:", same(a$fit_pred, b$fit_pred), "\n")
  if (!is.null(a$V) && !is.null(b$V) && !identical(a$V, b$V)) {
    dv <- sqrt(abs(diag(a$V))) - sqrt(abs(diag(b$V)))
    cat("   max |SE on - SE off|:", signif(max(abs(dv), na.rm = TRUE), 3),
        " NaN on/off:", sum(is.nan(diag(a$V))), sum(is.nan(diag(b$V))), "\n")
    print(rbind(on = sqrt(abs(diag(a$V))), off = sqrt(abs(diag(b$V)))))
  }
  if (!identical(a$re, b$re) && !is.null(a$re)) {
    cat("   max |ranef on - off|:", signif(max(abs(a$re - b$re)), 3), "\n")
  }
}

cat("\n== edge, seed 9031, fixed effects, se = FALSE\n")
a <- run("9031 fixed", y ~ x, dat(9031), TRUE)
b <- run("9031 fixed", y ~ x, dat(9031), FALSE)
cmp(a, b)
cat("\n== edge, seed 9031, fixed effects, se = TRUE\n")
a <- run("9031 fixed se", y ~ x, dat(9031), TRUE, se = TRUE)
b <- run("9031 fixed se", y ~ x, dat(9031), FALSE, se = TRUE)
cmp(a, b)
cat("\n== interior, seed 9001, fixed effects\n")
a <- run("9001 fixed", y ~ x, dat(9001), TRUE)
b <- run("9001 fixed", y ~ x, dat(9001), FALSE)
cmp(a, b)
for (s in c(9031, 9039, 9001)) {
  cat("\n== seed", s, "with (1 | g), 15 groups\n")
  d <- dat(s, G = 15)
  a <- run(paste(s, "(1|g)"), y ~ x + (1 | g), d, TRUE)
  b <- run(paste(s, "(1|g)"), y ~ x + (1 | g), d, FALSE)
  cmp(a, b)
}
cat("\n== weights and a predicted nu, seed 9031\n")
d <- dat(9031)
d$w <- rep(c(0.5, 1.5), 150)
a <- run("9031 weights", y | weights(w) ~ x, d, TRUE)
b <- run("9031 weights", y | weights(w) ~ x, d, FALSE)
cmp(a, b)
a <- run("9031 nu ~ x", bf(y ~ x, nu ~ x), d, TRUE)
cat("\n== right censored at the 90th percentile, seed 9031\n")
d <- dat(9031)
q9 <- quantile(d$y, 0.9)
d$cc <- ifelse(d$y > q9, "right", "none")
d$y <- pmin(d$y, q9)
a <- run("9031 cens", y | cens(cc) ~ x, d, TRUE)
b <- run("9031 cens", y | cens(cc) ~ x, d, FALSE)
cmp(a, b)
cat("\ndone\n")
