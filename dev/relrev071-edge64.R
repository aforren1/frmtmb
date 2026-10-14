# Release reviewer of 0.71.0: test-se-reparam.R's reparam_edge_fit(6, 4)
# ("a fit at the maximum with an sd at zero is not unconverged") on a
# build: code, warnings, lost set, logLik against lm() and lmer() (ML),
# what a larger sd gains, and the same on data moved by 0 to 4 ulps
# (seeds 1 to 20), to see whether the asserted outcome is decided by
# rounding.
#   Rscript dev/relrev071-edge64.R <lib> [nmoved]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "|", R.home(), "\n")
nm <- if (length(a) > 1) as.integer(a[2]) else 20L
set.seed(4)
d0 <- data.frame(x = stats::runif(300, 1, 2),
                 g = factor(rep(1:20, length.out = 300)))
d0$y <- sin(3 * d0$x) + stats::rnorm(300, 0, 0.3)
xn <- paste0("x", 1:6)
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "), "+ (1 | g)"))
one <- function(d, tag, verbose = FALSE) {
  for (k in 1:6) d[[xn[k]]] <- d$x^k
  w <- character()
  f <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  lmf <- stats::lm(stats::update(fo, . ~ . - (1 | g)), data = d)
  k <- length(stats::coef(lmf))
  ref <- unname(sqrt(diag(stats::vcov(lmf)) * (300 - k) / 300))
  se <- unname(sqrt(diag(suppressWarnings(vcov(f)))))[seq_len(k)]
  l <- ns$sdr_of(f)$se_lost
  gain <- tryCatch(ns$se_sd_gain_up(f, 1L), error = function(e) NA)
  cat(sprintf(paste0("%-6s code %d | dLL vs lm %+.2e | sd %.1e | lost %s |",
                     " gain up %.2e | max relSE %.2e | conv warn %d | grad",
                     " warn %d | other %d\n"),
              tag, f$opt$convergence,
              as.numeric(logLik(f)) - as.numeric(logLik(lmf)),
              exp(f$opt$par[["theta"]]),
              if (length(l)) paste(names(l), l, sep = ":", collapse = ",")
              else "-", as.numeric(gain), max(abs(se / ref - 1)),
              sum(grepl("did not report convergence", w)),
              sum(grepl("Large maximum absolute gradient", w)),
              sum(!grepl("did not report convergence|Large maximum", w))))
  if (verbose) for (x in w) cat("   W:", substr(gsub("\n", " ", x), 1, 300), "\n")
}
one(d0, "orig", verbose = TRUE)
m <- suppressMessages(lme4::lmer(stats::as.formula(paste("y ~",
       paste(sprintf("I(x^%d)", 1:6), collapse = " + "), "+ (1 | g)")),
       data = d0, REML = FALSE))
cat("lmer sd_g", as.data.frame(lme4::VarCorr(m))$sdcor[1], "\n")
for (s in seq_len(nm)) {
  set.seed(9100 + s)
  d <- d0
  for (v in c("x", "y")) d[[v]] <- d[[v]] * (1 + sample(0:4, 300, TRUE) * 2^-52)
  one(d, paste0("m", s))
}
