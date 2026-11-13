# The one boundary message dev/quadseam-setier-quad.R gave where
# glmer(nAGQ = 25) is not singular (arm A, binomial, sd 0): which seed,
# and what Laplace frmtmb, glmer Laplace and glmer nAGQ = 25 say there.
#   Rscript dev/quadseam-fa.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib: frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
gen <- function(s) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:60, each = 6)), x = rnorm(360))
  eta <- 0.2 + 0.5 * d$x + rnorm(60, 0, 0)[d$g]
  d$k <- rbinom(360, 10, plogis(eta)); d$nk <- 10 - d$k
  d
}
msgs <- function(expr) {
  m <- character()
  v <- withCallingHandlers(expr, message = function(x) {
    m <<- c(m, class(x)[1]); invokeRestart("muffleMessage")
  }, warning = function(x) {
    m <<- c(m, paste0("W:", class(x)[1])); invokeRestart("muffleWarning")
  })
  list(v = v, m = m)
}
for (s in 1:20) {
  d <- gen(s)
  q <- msgs(frm(cbind(k, nk) ~ x + (1 | g), family = binomial(), data = d,
                quadrature = TRUE))
  l25 <- suppressMessages(suppressWarnings(
    glmer(cbind(k, nk) ~ x + (1 | g), data = d, family = binomial,
          nAGQ = 25)))
  bq <- "frmtmb_boundary_fit" %in% q$m
  if (bq == isSingular(l25)) next
  lp <- msgs(frm(cbind(k, nk) ~ x + (1 | g), family = binomial(), data = d))
  l1 <- suppressMessages(suppressWarnings(
    glmer(cbind(k, nk) ~ x + (1 | g), data = d, family = binomial)))
  cat(sprintf(paste0("SEED %d quad: sd %.4g ll %.8f [%s] | laplace: sd %.4g",
                     " ll %.8f [%s]\n   glmer25 sd %.4g ll %.8f singular %s",
                     " | glmer1 sd %.4g ll %.8f singular %s\n"),
              s, exp(q$v$estimates$theta[1]), as.numeric(logLik(q$v)),
              paste(q$m, collapse = ","), exp(lp$v$estimates$theta[1]),
              as.numeric(logLik(lp$v)), paste(lp$m, collapse = ","),
              attr(VarCorr(l25)$g, "stddev"), as.numeric(logLik(l25)),
              isSingular(l25), attr(VarCorr(l1)$g, "stddev"),
              as.numeric(logLik(l1)), isSingular(l1)))
  # the quadrature fit's profile along its log sd, by rules taped there
  fn <- frmtmb:::quad_reanchored_fn(q$v)
  if (!is.null(fn)) {
    p <- q$v$opt$par
    for (sd in c(1e-6, 1e-4, 1e-3, 0.003, 0.01, 0.03, 0.1)) {
      pp <- p; pp[3] <- log(sd)
      cat(sprintf("   quad profile sd %.0e  -objective %.10f\n", sd,
                  -fn(pp)))
    }
  }
}
cat("END\n")
