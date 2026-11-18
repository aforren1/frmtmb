# Short review of lane alapse, fix round 1: a bound on the quantile's or
# sigma's coefficient (set_prior(lb =, ub =)) against the quantile
# search, on design B seed 109 (global p 0.7247 at log sigma -1.2668,
# logit 0.968; local p 0.4197 at log sigma -1.041, logit -0.324).
#   Rscript dev/reviews/alapse-short-bounds.R <lane|base>
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-alapse-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE, width = 150)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
cap <- function(expr) {
  wc <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(z) { wc <<- c(wc, class(z)[1])
      invokeRestart("muffleWarning") },
    message = function(z) invokeRestart("muffleMessage"))
  list(v = v, wc = wc)
}
set.seed(2000 + 109)
dB <- data.frame(x = rnorm(100))
dB$y <- 1 + 0.5 * dB$x + rnorm(100)
cases <- list(
  "quantile Int in [-2, 0.3]" = list(
    fo = bf(y ~ x, quantile ~ 1),
    pr = set_prior("", class = "Intercept", dpar = "quantile",
                   lb = -2, ub = 0.3)),
  "sigma Int in [-1.1, 2]" = list(
    fo = bf(y ~ x, sigma ~ 1),
    pr = set_prior("", class = "Intercept", dpar = "sigma",
                   lb = -1.1, ub = 2)),
  "quantile Int in [-2, 0.3], fixed-start" = list(
    fo = bf(y ~ x, quantile ~ 1),
    pr = set_prior("", class = "Intercept", dpar = "quantile",
                   lb = -2, ub = 0.3), ctl = frmtmb_control(restarts = 0)))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  r <- cap(frm(cs$fo, family = asym_laplace(), data = dB, prior = cs$pr,
               control = cs$ctl %||% frmtmb_control()))
  if (inherits(r$v, "error")) {
    cat(nm, "ERROR", conditionMessage(r$v), "\n"); next
  }
  f <- r$v
  bx <- f$cache$bounds
  np <- length(f$opt$par)
  lo <- rep_len(bx$lower %||% -Inf, np)
  hi <- rep_len(bx$upper %||% Inf, np)
  out <- f$opt$par < lo | f$opt$par > hi
  dp <- ns$eval_dpars(f)[[1]]
  cat(sprintf(paste0("%-40s code %d global %-5s | p %.4f sigma %.4f logLik ",
                     "%.6f | par %s | box [%s] .. [%s] | OUTSIDE BOX: %s | ",
                     "warn %s\n"),
              nm, f$opt$convergence, isTRUE(f$opt$ald_cert$global),
              dp$quantile[1], dp$sigma[1], as.numeric(logLik(f)),
              paste(sprintf("%.4f", f$opt$par), collapse = " "),
              paste(format(lo, digits = 3), collapse = " "),
              paste(format(hi, digits = 3), collapse = " "),
              paste(names(f$opt$par)[out], collapse = ",") ,
              paste(r$wc, collapse = ",")))
}
cat("DONE\n")
