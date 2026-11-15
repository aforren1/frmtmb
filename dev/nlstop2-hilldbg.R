# Lane nlstop2, item 5093: the Hill body on step data (the review's
# fndebug row): the objective gain along the flattest direction in log
# magnitudes at each factor, before and after one Newton step across.
# Usage: Rscript dev/nlstop2-hilldbg.R
source("dev/nlstop2-common.R")
ns <- asNamespace("frmtmb")
set.seed(1)
n <- 1000
d <- data.frame(x = c(runif(n / 2, 0.5, 4.5), runif(n / 2, 5.5, 10)))
d$y <- 10 * (d$x > 5) + rnorm(n, 0, 0.5)
r <- tt(frm(bf(y ~ Emax * x^h / (EC50^h + x^h), Emax + EC50 + h ~ 1,
               nl = TRUE), data = d, start = list(beta = c(10, 5, 2))))
cat("fit", llv(r), r$opt$convergence, wsum(r), "\n")
obj <- r$obj; p0 <- r$opt$par; j <- which(names(p0) == "beta")
g0 <- drop(obj$gr(p0))[j]
H <- matrix(NA, 3, 3)
for (i in 1:3) {
  h <- 1e-4 * abs(p0[j[i]]); up <- p0; dn <- p0
  up[j[i]] <- up[j[i]] + h; dn[j[i]] <- dn[j[i]] - h
  H[, i] <- (drop(obj$gr(up))[j] - drop(obj$gr(dn))[j]) / (2 * h)
}
D <- p0[j]
Hu <- (H + t(H)) / 2 * outer(D, D) + diag(D * g0, 3)
e <- eigen(Hu, symmetric = TRUE)
v <- e$vectors[, 3]
for (fac in c(1000, 100, 10, 4, 2, 1.5)) {
  for (s in c(1, -1)) {
    t <- log(fac) / max(abs(v))
    q <- p0; q[j] <- p0[j] * exp(s * t * v)
    f1 <- obj$fn(q)
    gq <- drop(obj$gr(q))[j]
    E <- e$vectors[, 1:2]
    du <- -drop(E %*% (crossprod(E, q[j] * gq) / e$values[1:2]))
    q2 <- q; q2[j] <- q[j] * exp(du)
    f2 <- obj$fn(q2)
    cat(sprintf("fac %g side %d gate %.3g gain %.6g after newton %.6g est %s\n",
                fac, s, 0.5 * e$values[3] * t^2, r$opt$objective - f1,
                r$opt$objective - f2,
                paste(format(q2[j], digits = 5), collapse = ",")))
  }
}
for (fac in c(1000, 100, 10, 4, 2, 1.5)) {
  for (s in c(1, -1)) {
    q <- p0; q[j[3]] <- p0[j[3]] * fac^s
    cat(sprintf("h alone x %g^%d: gain %.6g\n", fac, s,
                r$opt$objective - obj$fn(q)))
  }
}
# the profile in h: every other outer parameter optimized with h held
for (fac in c(1.5, 2, 4, 8)) {
  hv <- p0[j[3]] * fac
  fr <- setdiff(seq_along(p0), j[3])
  pf <- function(z) { q <- p0; q[fr] <- z; q[j[3]] <- hv; obj$fn(q) }
  o <- nlminb(p0[fr], pf)
  cat(sprintf("profile h x %g: gain %.6g code %d\n", fac,
              r$opt$objective - o$objective, o$convergence))
}
cat("nl_limit_message():", ns$nl_limit_message(r), "\n")
