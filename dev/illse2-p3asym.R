# Lane illse2, item 3: the asymmetry of the Hessian in the design's (or
# Jacobian's) coordinates, max |H - H'| / 2 over max |H|, on the nl
# cases of dev/illse2-p3nl.R and on gaussian raw polynomials, by step.
#   ILLSE2_ARM=lane Rscript dev/illse2-p3asym.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
mk <- list(
  rat = function(seed) {
    set.seed(seed)
    d <- data.frame(x = stats::runif(200, 1, 2))
    d$y <- sin(3 * d$x) / (1 + 0.3 * d$x) + stats::rnorm(200, 0, 0.3)
    list(d = d, f = bf(y ~ (b0 + b1 * x + b2 * x^2 + b3 * x^3 + b4 * x^4 +
                              b5 * x^5 + b6 * x^6) / (1 + dd * x),
                       b0 + b1 + b2 + b3 + b4 + b5 + b6 + dd ~ 1, nl = TRUE))
  },
  sin = function(seed) {
    set.seed(seed)
    d <- data.frame(x = stats::runif(300, 1, 2), z = stats::runif(300, 0, 10))
    d$y <- sin(3 * d$x) + 1.5 * sin(1.3 * d$z) + stats::rnorm(300, 0, 0.3)
    list(d = d, f = bf(y ~ b0 + b1 * x + b2 * x^2 + b3 * x^3 + b4 * x^4 +
                         b5 * x^5 + b6 * x^6 + A * sin(w * z),
                       b0 + b1 + b2 + b3 + b4 + b5 + b6 + A + w ~ 1, nl = TRUE),
         start = list(beta = c(unname(stats::coef(stats::lm(
           d$y ~ outer(d$x, 1:6, `^`)))), 1, 1.3)))
  },
  poly8 = function(seed) {
    set.seed(seed)
    d <- data.frame(x = stats::runif(200, 1, 2))
    d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
    for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
    list(d = d, f = bf(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8))
  })
for (cs in names(mk)) for (seed in 1:3) {
  s <- mk[[cs]](seed)
  fit <- suppressWarnings(if (is.null(s$start)) frm(s$f, data = s$d) else
    frm(s$f, data = s$d, start = s$start))
  p <- fit$opt$par
  M <- ns$reparam_basis(fit$frame, names(p), est = fit$estimates)
  if (is.null(M)) { cat(cs, seed, "no basis\n"); next }
  for (h in c(1e-3, 1e-4)) {
    hh <- ns$reparam_hessian(fit$obj, p, M, step = h)
    cat(sprintf("%-5s seed %d step %g | asym %.2e | min diag %.2e\n", cs,
                seed, h, max(hh$E) / max(abs(hh$H)), min(diag(hh$H))))
  }
}
