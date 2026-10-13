# Lane seflags, punch round 2b: degree 8 without random effects, seed 7.
# Is the 4 percent excess of the SEs over the exact ones (at the fit's
# sigma) the coupling of the "illcond" direction with log sigma at a
# point the optimizer left short along that direction? A Newton step
# along the direction, then the same covariance (span of the two
# smallest re-measured) and the cross term, before and after.
#   Rscript dev/seflags-p2b-span2.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
set.seed(7)
d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] + stats::rnorm(200, 0, 0.2)
for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
rhs <- paste0("x", 1:8, collapse = " + ")
f <- suppressWarnings(frm(stats::as.formula(paste("y ~", rhs)), data = d))
on <- ns$outer_par_names(f)
X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
j <- match(colnames(X), on)
an <- function(p) {
  H <- f$obj$he(p); H <- (H + t(H)) / 2
  free <- seq_along(p); Df <- sqrt(abs(diag(H)))
  e <- eigen(H / outer(Df, Df), symmetric = TRUE)
  n <- length(e$values); S <- (n - 1):n
  sp <- ns$se_span_remeasure(f, p, free, e, Df, S)
  U <- e$vectors; ev <- e$values
  U[, S] <- sp$vectors; ev[S] <- sp$values
  V <- U %*% (t(U) / ev) / outer(Df, Df)
  sg <- exp(p[grepl("^sigma", on)])
  q <- qr(X / sg, tol = 1e-30)
  Ri <- backsolve(qr.R(q), diag(ncol(X)))
  ex <- numeric(ncol(X)); ex[q$pivot] <- sqrt(rowSums(Ri^2))
  dn <- sp$vectors[, which.min(sp$values)] / Df
  g <- as.numeric(f$obj$gr(p))
  ks <- grep("^sigma", on)
  cat(sprintf("  grad . d %.3g | H[d, log sigma] %.3g | SE/exact %s\n",
              sum(g * dn), sum(H[ks, ] * dn),
              paste(range(signif(sqrt(diag(V))[j] / ex, 4)), collapse = " to ")))
  list(dn = dn, c = min(sp$values), g = g)
}
p <- f$opt$par
cat("at the fit\n"); a <- an(p)
for (it in 1:3) {
  p <- p - sum(a$g * a$dn) / a$c * a$dn
  cat("after Newton step", it, "along d, objective", f$obj$fn(p) - f$opt$objective, "\n")
  a <- an(p)
}
