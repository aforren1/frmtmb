# Release reviewer of 0.71.0: does the rewritten start of
# test-se-conv.R's last test still pin lane illse's R2-1 fix (the
# headroom asks the likelihood along an empty sd row)? The test's fit on
# the given library, as installed and with se_sd_gain_up() stubbed to 0
# (the headroom without the fix), printing what the test asserts on.
#   Rscript dev/relrev071-seconv.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
conv_conds <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}
one <- function(tag, start) {
  set.seed(1)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:20, length.out = n)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.5)[d$g] +
    stats::rnorm(n, 0, 0.4)
  r <- conv_conds(frm(y ~ x + (1 | g), data = d, start = start,
                      control = frmtmb_control(optCtrl = list(
                        iter.max = 3, eval.max = 1000))))
  f <- r$value
  m <- suppressMessages(lme4::lmer(y ~ x + (1 | g), data = d, REML = FALSE))
  cat(sprintf("%-28s code %d | theta %.3f | gap to lmer %.4f | sd grad %.2e\n",
              tag, f$opt$convergence, f$opt$par[["theta"]],
              as.numeric(logLik(f)) - as.numeric(logLik(m)),
              drop(f$obj$gr(f$opt$par))[length(f$opt$par)]))
  for (w in r$warnings) cat("   W:", substr(gsub("\n", " ", w), 1, 230), "\n")
}
st_new <- list(beta = c(0, 0), betad = 0, theta = -15)
st_old <- list(theta = -15)
one("new start, as installed", st_new)
one("old start, as installed", st_old)
if (exists("se_sd_gain_up", ns)) {
  orig <- get("se_sd_gain_up", ns)
  unlockBinding("se_sd_gain_up", ns)
  assign("se_sd_gain_up", function(fit, k, ...) 0, envir = ns)
  one("new start, gain stubbed 0", st_new)
  one("old start, gain stubbed 0", st_old)
  assign("se_sd_gain_up", orig, envir = ns)
}
# a start at lm()'s maximum given sd = 0, no iteration: only the sd row
# is short, which is R2-1's premise whatever start frm() would choose
set.seed(1)
n <- 300
d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:20, length.out = n)))
d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.5)[d$g] + stats::rnorm(n, 0, 0.4)
lm0 <- stats::lm(y ~ x, data = d)
st_lm <- list(beta = unname(coef(lm0)),
              betad = log(sqrt(mean(residuals(lm0)^2))), theta = -15)
for (it in c(0, 1)) {
  ctl <- frmtmb_control(optCtrl = list(iter.max = it, eval.max = 1000))
  r <- conv_conds(frm(y ~ x + (1 | g), data = d, start = st_lm, control = ctl))
  cat(sprintf("lm start iter.max %d, as installed: code %d theta %.2f\n", it,
              r$value$opt$convergence, r$value$opt$par[["theta"]]))
  for (w in r$warnings) cat("   W:", substr(gsub("\n", " ", w), 1, 200), "\n")
  if (exists("se_sd_gain_up", ns)) {
    assign("se_sd_gain_up", function(fit, k, ...) 0, envir = ns)
    r <- conv_conds(frm(y ~ x + (1 | g), data = d, start = st_lm, control = ctl))
    cat(sprintf("lm start iter.max %d, gain stubbed: code %d\n", it,
                r$value$opt$convergence))
    for (w in r$warnings) cat("   W:", substr(gsub("\n", " ", w), 1, 200), "\n")
    assign("se_sd_gain_up", orig, envir = ns)
  }
}
