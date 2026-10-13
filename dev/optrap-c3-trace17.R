# Lane optrap, consolidation fixes: the Newton steps of degree 7, seed
# 17, no group term (dev/optrap-c3-repoly.R's data), traced: at each
# step the predicted and actual gain, whether the Hessian factors, the
# natural-unit gradient. Which rule stops the steps short of the
# maximum.
#   Rscript dev/optrap-c3-trace17.R <lib> [deg] [seed]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
deg <- if (length(a) > 1) as.integer(a[2]) else 7L
s <- if (length(a) > 2) as.integer(a[3]) else 17L
set.seed(1000 * deg + s)
d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + ")))
ns <- asNamespace("frmtmb")
utils::assignInNamespace("newton_polish", function(obj, opt, ...) opt,
                         "frmtmb")
f <- suppressWarnings(frm(fo, data = d))
X <- cbind(1, as.matrix(d[paste0("x", 1:deg)]))
r <- qr.qty(qr(X, LAPACK = TRUE), d$y)[-seq_len(ncol(X))]
best <- -0.5 * 200 * (log(2 * pi) + log(sum(r^2) / 200) + 1)
p <- f$opt$par
u <- f$par_units %||% rep(1, length(p))
v <- f$opt$objective
cat("nlminb stop: short", -v - best, "\n")
for (it in 1:8) {
  g <- drop(f$obj$gr(p)) * u
  H <- f$obj$he(p) * outer(u, u)
  R <- tryCatch(chol((H + t(H)) / 2), error = function(e) NULL)
  ev <- range(eigen((H + t(H)) / 2, only.values = TRUE)$values)
  if (is.null(R)) { cat("step", it, "no Cholesky; eigen range", ev, "\n"); break }
  st <- -drop(backsolve(R, forwardsolve(t(R), g)))
  gain <- -sum(g * st) / 2
  q <- p + st * u
  vq <- f$obj$fn(q)
  cat(sprintf("step %d | |g| %.3g | eigen %.3g..%.3g | predicted %.4g actual %.4g ratio %.3f | short after %.4g\n",
              it, max(abs(g)), ev[1], ev[2], gain, v - vq, (v - vq) / gain,
              -vq - best))
  if (!(vq < v)) break
  p <- q; v <- vq
}
