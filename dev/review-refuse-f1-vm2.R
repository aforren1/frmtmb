# Short check, item 6: the von_mises line against brms's window when a reader
# reads the value linearly. Same data as the laplace script's first cell (seed
# 11200 + 7), then 12 replicates. Three fits: frm() on the line; frm() with the
# window chart restored by a mutant of mi_ml_frame() in this process; the exact
# maximum on the window. Which row's value leaves the window?
source("dev/review-refuse-lib.R")
rvm <- function(mu, k) {
  vapply(mu, function(m) { repeat { w <- runif(1, -pi, pi); if (runif(1) < exp(k * (cos(w) - 1))) break }; ((w + m + pi) %% (2 * pi)) - pi }, 0)
}
exact <- function(d, start) {
  rows <- which(is.na(d$y)); ob <- setdiff(seq_len(nrow(d)), rows)
  vals <- seq(-pi, pi, length.out = 1001); dg <- diff(vals)[1]
  ldv <- function(y, eta, lk) exp(lk) * cos(y - 2 * atan(eta)) - log(2 * pi) - log(besselI(exp(lk), 0))
  nll <- function(p) {
    ey <- p[1] + p[2] * d$x
    s <- sum(ldv(d$y[ob], ey[ob], p[3])) + sum(dnorm(d$z[ob], p[4] + p[5] * d$y[ob] + p[6] * d$x[ob], exp(p[7]), log = TRUE))
    for (i in rows) {
      li <- ldv(vals, ey[i], p[3]) + dnorm(d$z[i], p[4] + p[5] * vals + p[6] * d$x[i], exp(p[7]), log = TRUE)
      mx <- max(li); s <- s + mx + log(sum(exp(li - mx)) * dg)
    }
    if (!is.finite(s)) 1e10 else -s
  }
  o <- nlminb(start, nll, control = list(iter.max = 2000, eval.max = 4000, rel.tol = 1e-10))
  list(par = o$par, ll = -o$objective, se = sqrt(diag(solve(stats::optimHess(o$par, nll)))), fn = nll)
}
orig <- get("mi_ml_frame", asNamespace("frmtmb"))
window_mi_ml <- orig
body(window_mi_ml) <- do.call(substitute, list(body(orig), list()))
txt <- deparse(orig)
i <- grep('identical(m[["chart"]], "circle")', txt, fixed = TRUE)
stopifnot(length(i) == 1L)
txt[i] <- gsub('identical(m[["chart"]], "circle")', "FALSE", txt[i], fixed = TRUE)
window_mi_ml <- eval(parse(text = txt)); environment(window_mi_ml) <- asNamespace("frmtmb")
fitit <- function(d, line = TRUE) {
  assignInNamespace("mi_ml_frame", if (line) orig else window_mi_ml, "frmtmb")
  on.exit(assignInNamespace("mi_ml_frame", orig, "frmtmb"))
  tr(frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y) + x) + gaussian(), data = d))
}
pars <- function(a) {
  e <- a$estimates; c(e$beta[grep("^y_", names(e$beta))], e$betad[grep("^y_", names(e$betad))],
                      e$beta[grep("z_.Intercept", names(e$beta))], e$beta[grep("miy", names(e$beta))], e$beta[grep("^z_x", names(e$beta))], e$betad[grep("^z_", names(e$betad))])
}
set.seed(11200 + 7)
n <- 200
res <- list()
for (r in 0:12) {
  x <- rnorm(n); eta <- 0.3 + 0.5 * x
  y <- rvm(2 * atan(eta), 3)
  d <- data.frame(x = x, y = y); d$z <- 0.5 + 0.7 * y + 0.3 * x + rnorm(n, 0, 0.7)
  rows <- sample(n, 60); ytrue <- d$y[rows]; d$y[rows] <- NA
  a <- fitit(d, TRUE); w <- fitit(d, FALSE)
  if (is_err(a) || is_err(w)) { pr("r=%d line %s window %s", r, msg(a, 60), msg(w, 60)); next }
  pa <- pars(a); pw <- pars(w)
  ex <- tr(exact(d, pw))
  if (is_err(ex)) { pr("r=%d exact %s", r, msg(ex, 60)); next }
  ex2 <- tr(exact(d, pa))
  if (!is_err(ex2) && ex2$ll > ex$ll) ex <- ex2
  v <- frmtmb:::mi_values(a, "y")[sort(rows)]
  sh <- function(p) (p - ex$par) / ex$se
  pr("r=%2d line: leaves the window %2d of 60, range %.2f %.2f; max|shift| line %.2f (zy %+.2f) window %.2f (zy %+.2f); logLik line %.2f window %.2f exact %.2f; zy line %.3f window %.3f exact %.3f",
     r, sum(abs(v) >= pi), min(v), max(v), max(abs(sh(pa))), sh(pa)[5], max(abs(sh(pw))), sh(pw)[5], as.numeric(logLik(a)), as.numeric(logLik(w)), ex$ll, pa[5], pw[5], ex$par[5])
}
