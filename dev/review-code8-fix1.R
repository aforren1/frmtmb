# Review of lane code8, item 3: brms fixture 1. Are the sigma SEs the
# lane reports (0.2848, 0.4831) those of the boundary maximum? Reference:
# a Richardson Hessian of the Laplace objective over every parameter but
# the correlation (theta_3) and the smooth's sd (theta_4), pseudo-inverted
# past the exact alias, at the stop and at points with the correlation
# parameter held at other values and the rest re-optimized.
#   RV_ARM=after Rscript dev/review-code8-fix1.R
source("dev/review-code8-common.R")
options(width = 180)
fit0 <- get("brmsfit_example1", envir = asNamespace("brms"))
d <- as.data.frame(fit0$data)
attr(d, "terms") <- NULL
form <- bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) +
             (1 + Trt | visit) + arma(visit, patient), sigma ~ Trt)
r <- rfit(form, family = student(), data = d)
f <- r$f
cat("conds:", cond_txt(r), "\n")
p <- f$opt$par
nm <- ns$outer_par_names(f)
s <- se_of(f)
print(s$lost)
isg <- match(c("sigma_(Intercept)", "sigma_Trt1"), nm)
cat("lane SEs sigma:", s$se[isg], "\n")
i3 <- match("theta_3", nm); i4 <- match("theta_4", nm)
cat("theta_3", p[i3], "theta_4", p[i4], "\n")
obj <- f$obj
keep <- setdiff(seq_along(p), c(i3, i4))
refse <- function(P) {
  gk <- function(q) { z <- P; z[keep] <- q; drop(obj$gr(z))[keep] }
  J <- numDeriv::jacobian(gk, P[keep], method = "Richardson")
  H <- (J + t(J)) / 2
  D <- sqrt(abs(diag(H)))
  e <- eigen(H / outer(D, D), symmetric = TRUE)
  ok <- e$values > 1e-8 * max(e$values)
  Vu <- e$vectors[, ok] %*% (t(e$vectors[, ok]) / e$values[ok])
  V <- Vu / outer(D, D)
  se <- sqrt(diag(V))
  names(se) <- nm[keep]
  list(se = se, ev = e$values, ndrop = sum(!ok))
}
r0 <- refse(p)
cat("at stop: dropped", r0$ndrop, "smallest ev", signif(tail(r0$ev, 3), 3), "\n")
cat("ref SEs sigma:", r0$se[c("sigma_(Intercept)", "sigma_Trt1")], "\n")
cmp <- data.frame(par = nm[keep], lane = signif(s$se[keep], 5),
                  ref = signif(r0$se, 5))
print(cmp)
for (t3 in c(-30, -100, -300, -3000)) {
  st <- p; st[i3] <- t3
  oth <- setdiff(seq_along(p), c(i3, i4))
  fo <- function(q) { z <- st; z[oth] <- q; obj$fn(z) }
  go <- function(q) { z <- st; z[oth] <- q; drop(obj$gr(z))[oth] }
  o <- nlminb(p[oth], fo, go, control = list(iter.max = 500, eval.max = 1000,
                                              rel.tol = 1e-14))
  P <- st; P[oth] <- o$par
  rr <- refse(P)
  cat(sprintf("theta_3 %6g: objective %.8f  sigma SEs %.5f %.5f  dropped %d\n",
              t3, o$objective, rr$se[["sigma_(Intercept)"]],
              rr$se[["sigma_Trt1"]], rr$ndrop))
}
cat("DONE\n")
