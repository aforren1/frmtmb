# Lane boundfam, item 3098: three-dimensional us() blocks whose covariance
# is singular (rank 2), so a Cholesky row runs off with no correlation at
# +/-1. Which seeds give a zero-width correlation interval (se_t == 0)?
#   Rscript dev/boundfam-us3.R <lib or "base"> <seed from> <seed to>
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
seeds <- seq(as.integer(a[2]), as.integer(a[3]))
for (s in seeds) {
  set.seed(s)
  ng <- 30
  m <- 10
  d <- data.frame(x = rnorm(ng * m), z = rnorm(ng * m),
                  g = factor(rep(seq_len(ng), each = m)))
  u0 <- rnorm(ng, 0, 0.8)
  u1 <- rnorm(ng, 0, 0.5)
  # the z slope is a combination of the other two: rank 2, and no
  # correlation at +/-1
  u2 <- 0.4 * u0 - 0.6 * u1
  d$y <- 1 + 0.5 * d$x - 0.3 * d$z + u0[d$g] + u1[d$g] * d$x +
    u2[d$g] * d$z + rnorm(ng * m, 0, 0.5)
  w <- character(0)
  f <- tryCatch(withCallingHandlers(
    frm(y ~ x + z + (1 + x + z | g), data = d),
    warning = function(e) {
      w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning")
    },
    message = function(e) {
      w <<- c(w, paste("M:", conditionMessage(e)))
      invokeRestart("muffleMessage")
    }), error = function(e) NULL)
  if (is.null(f)) {
    cat(sprintf("seed %d: error\n", s))
    next
  }
  th <- f$estimates$theta
  lost <- tryCatch(ns$sdr_of(f)$se_lost, error = function(e) NULL)
  tr <- suppressWarnings(ns$varcorr_trans_rows(f))
  cv <- suppressWarnings(confint_varcorr(f))
  zw <- which(is.finite(cv$lwr) & cv$lwr == cv$upr)
  cat(sprintf("seed %d: code %d logLik %.6f theta %s lost %s\n", s,
              f$opt$convergence, as.numeric(logLik(f)),
              paste(signif(th, 5), collapse = " "),
              paste(names(lost), lost, collapse = ", ")))
  cat("   se_t:", paste(signif(tr$se_t, 4), collapse = " "), "\n")
  cat("   zero-width rows:", if (length(zw)) {
    paste(cv$term[zw], signif(cv$estimate[zw], 6), collapse = "; ")
  } else "none", "\n")
  for (x in w) cat("   W:", substr(gsub("\\s+", " ", x), 1, 300), "\n")
}
cat("done\n")
