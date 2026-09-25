# Punch round 1, M2: is the link-scale se.fit a detector for an absurd
# predict() row, once legitimate extrapolation is in the calibration
# set? Recorded per row: se_eta, and the row's leverage in the
# training design, h = n * x0' (X'X)^-1 x0 (1 on average in sample),
# which says whether the row is extrapolation in design space.
#   PREDFIX_ARM=base Rscript dev/predfix-p1-se.R > dev/predfix-log/p1-se.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
q <- function(expr) suppressWarnings(suppressMessages(expr))
lev <- function(f, nd) {
  X <- as.matrix(f$frame$linpreds[[1]]$X)
  tt <- stats::delete.response(f$frame$linpreds[[1]]$terms)
  X0 <- stats::model.matrix(tt, stats::model.frame(tt, nd,
                             xlev = f$frame$linpreds[[1]]$xlevels))
  XtXi <- tryCatch(solve(crossprod(X)), error = function(e) MASS::ginv(crossprod(X)))
  nrow(X) * rowSums((X0 %*% XtXi) * X0)
}
rows <- list()
add <- function(label, kind, f, nd) {
  se <- q(frm_linpred(f, newdata = nd, se.fit = TRUE))$se.fit
  rows[[length(rows) + 1L]] <<- data.frame(label = label, kind = kind,
    row = seq_len(nrow(nd)), se_eta = se, leverage = lev(f, nd))
}
for (seed in 1:10) {
  set.seed(seed)
  n <- 200
  x <- rnorm(n)
  for (s in c(4.5, 5.5, 6.5)) {
    d <- data.frame(x = x, y = exp(rnorm(n, 1 + 0.3 * x, s)))
    add(paste("lognormal sigma", s), "genuine",
        q(frm(bf(y ~ x) + lognormal(), data = d)), data.frame(x = c(0, 2)))
  }
  d <- data.frame(x = x, y = rnbinom(n, mu = exp(2 + 0.3 * x), size = 0.05))
  add("negbinomial shape 0.05", "genuine",
      q(frm(bf(y ~ x) + negbinomial(), data = d)), data.frame(x = c(0, 2)))
  d <- data.frame(x = x, y = 1 + 0.3 * x + rt(n, 1.05))
  add("student t(1.05)", "genuine", q(frm(bf(y ~ x) + student(), data = d)),
      data.frame(x = c(0, 2)))
  d <- data.frame(x = x, y = rpois(n, exp(0.3 + 0.4 * x)))
  add("poisson extrapolation", "extrapolation",
      q(frm(bf(y ~ x) + poisson(), data = d)), data.frame(x = c(30, 300)))
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  yy <- rpois(n, exp(0.5 + 0.3 * x))
  yy[g == "c"] <- 0L
  add("poisson all-zero cell", "flat",
      q(frm(bf(y ~ x + g) + poisson(), data = data.frame(x = x, g = g, y = yy))),
      data.frame(x = 0, g = factor(c("a", "c"), levels = c("a", "b", "c"))))
  for (e in c(1e-3, 1e-2)) {
    x2 <- x + rnorm(n, 0, e)
    add(paste("poisson x2 = x + N(0,", e, ")"), "flat",
        q(frm(bf(y ~ x + x2) + poisson(),
              data = data.frame(x = x, x2 = x2,
                                y = rpois(n, exp(0.5 + 0.3 * x))))),
        data.frame(x = c(0, 1), x2 = c(0, -1)))
  }
}
r <- do.call(rbind, rows)
saveRDS(r, "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/p1-se.rds")
options(width = 200)
key <- paste(r$label, "row", r$row)
tab <- data.frame(kind = tapply(r$kind, key, `[`, 1),
                  se_min = tapply(r$se_eta, key, min),
                  se_max = tapply(r$se_eta, key, max),
                  lev_min = tapply(r$leverage, key, min),
                  lev_max = tapply(r$leverage, key, max))
print(tab[order(tab$kind, tab$se_max), ], digits = 3)
