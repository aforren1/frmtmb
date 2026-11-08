# Lane bigx, short check M2: products and powers of shifted columns, and
# raw polynomials far from zero, against the centered model (or glm() on
# the same columns): which columns the frame keeps, the condition it
# reports, and the fit.
#   BIGX_ARM=before|after Rscript dev/bigx-prod.R
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
conds <- function(expr) {
  w <- character(0)
  m <- character(0)
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(x) {
      w <<- c(w, paste0("[", class(x)[1], "] ",
                        substr(conditionMessage(x), 1, 110)))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, substr(conditionMessage(x), 1, 110))
      invokeRestart("muffleMessage")
    })
  list(v = v, w = w, m = m)
}
show <- function(lbl, r, llref) {
  if (inherits(r$v, "error")) {
    cat(lbl, "ERROR", conditionMessage(r$v), "\n")
    return(invisible())
  }
  f <- r$v
  cat(sprintf("%-40s kept %d code %d ll - ref %+.3e\n", lbl,
              length(f$estimates[["beta"]]), f$opt$convergence,
              as.numeric(logLik(f)) - llref))
  for (x in r$w) cat("    warning:", x, "\n")
  for (x in r$m) cat("    message:", x, "\n")
}
for (m in c(1e4, 1e5, 1e6)) for (s in 1:3) {
  set.seed(s)
  n <- 200
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  d <- data.frame(z1, z2, w1 = z1 + m, w2 = z2 + m)
  d$y <- 1 + 0.5 * z1 - 0.3 * z2 + 0.4 * z1 * z2 + rnorm(n, 0, 0.5)
  d$yp <- rpois(n, exp(0.3 + 0.4 * z1 - 0.2 * z1^2))
  rc <- frm(bf(y ~ z1 * z2), family = gaussian(), data = d)
  show(sprintf("gaussian w1 * w2, m %g, seed %d", m, s),
       conds(frm(bf(y ~ w1 * w2), family = gaussian(), data = d)),
       as.numeric(logLik(rc)))
  pc <- frm(bf(yp ~ z1 + I(z1^2)), family = poisson(), data = d)
  show(sprintf("poisson w1 + I(w1^2), m %g, seed %d", m, s),
       conds(frm(bf(yp ~ w1 + I(w1^2)), family = poisson(), data = d)),
       as.numeric(logLik(pc)))
}
for (s in 501:503) for (deg in c(4, 6, 8)) {
  set.seed(s)
  n <- 200
  xx <- runif(n, 1000, 1010)
  d <- data.frame(xx = xx, yy = sin(xx / 2) + rnorm(n, 0, 0.3))
  fo <- as.formula(sprintf("yy ~ poly(xx, %d, raw = TRUE)", deg))
  r <- conds(frm(bf(fo), family = gaussian(), data = d))
  llref <- if (inherits(r$v, "error")) NA else {
    X <- r$v$frame$linpreds[[1]]$X
    g <- glm.fit(X, d$yy)
    -n / 2 * (log(2 * pi * sum(g$residuals^2) / n) + 1)
  }
  show(sprintf("raw poly %d on [1000, 1010], seed %d (ref: glm, kept)",
               deg, s), r, llref)
  if (!inherits(r$v, "error")) {
    se <- suppressWarnings(sqrt(diag(vcov(r$v))))
    cat("    SEs finite", sum(is.finite(se)), "of", length(se), "\n")
  }
}
