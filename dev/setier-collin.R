# Lane setier: does the eigenvalue gate before tier 1 take standard
# errors from ill-conditioned but identified designs the field fits?
# Raw polynomials of growing degree on x in [1, 2] (no random effects:
# the exact Hessian), and with a random intercept (finite differences),
# against lm()'s and lme4's standard errors.
#   Rscript dev/setier-collin.R <lib or "base">
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
if (identical(args[1], "src")) {
  suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-setier",
                                     quiet = TRUE))
} else suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
for (deg in 3:7) {
  set.seed(deg)
  d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  fo <- stats::as.formula(paste0("y ~ ", paste0("I(x^", 1:deg, ")",
                                               collapse = " + ")))
  ref <- sqrt(diag(vcov(lm(fo, data = d))))
  X <- model.matrix(fo, d)
  S <- crossprod(X); S <- S / sqrt(outer(diag(S), diag(S)))
  ev <- eigen(S, only.values = TRUE)$values
  w <- character()
  f <- withCallingHandlers(frm(fo, data = d), warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")})
  se <- suppressWarnings(fixef(f)[, "Est.Error"])
  cat(sprintf(paste0("deg %d no RE: X'X unit-diag ev ratio %.2e; ",
                     "frmtmb/lm SE ratio range %s; lost %d; warnings %d\n"),
              deg, min(ev) / max(ev),
              paste(format(range(se / ref), digits = 4), collapse = " "),
              length(frmtmb:::sdr_of(f)$se_lost), length(w)))
  fo2 <- stats::update(fo, . ~ . + (1 | g))
  m <- lme4::lmer(fo2, data = d, REML = FALSE)
  ref2 <- sqrt(diag(as.matrix(vcov(m))))
  w <- character()
  f2 <- withCallingHandlers(frm(fo2, data = d), warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")},
    message = function(c) invokeRestart("muffleMessage"))
  se2 <- suppressWarnings(fixef(f2)[, "Est.Error"])
  cat(sprintf(paste0("deg %d (1 | g): frmtmb/lme4 SE ratio range %s; ",
                     "lost %d; warnings %d\n"), deg,
              paste(format(range(se2 / ref2), digits = 4), collapse = " "),
              length(frmtmb:::sdr_of(f2)$se_lost), length(w)))
}
