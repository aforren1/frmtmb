# Reviewer of lane setier: which directions se_tier3() probes on c0k seed
# 77 and what se_curvature_real() returns (traced).
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(77)
G <- 8
dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(G, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
trace("se_curvature_real", where = ns, print = FALSE,
      exit = quote(cat("se_curvature_real lambda", signif(lambda, 3),
                       "->", returnValue(), "\n")))
trace("se_line_probe", where = ns, print = FALSE,
      exit = quote(cat("se_line_probe lambda", signif(lambda, 3),
                       "->", returnValue(), "\n")))
f <- suppressMessages(suppressWarnings(frm(bf(yn ~ c0 + exp(a)^k,
  c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE), data = dn)))
print(ns$sdr_of(f)$se_lost)
an <- f$cache$se_analysis
print(an$ev)
# the probed direction, its losses, and the objective along the exact ridge
H <- f$cache$hessian_fixed$H
nm <- ns$outer_par_names(f)
D <- sqrt(abs(diag(H)))
e <- eigen(H / outer(D, D), symmetric = TRUE)
k <- 5
p <- f$opt$par
tol <- 1e-3
for (mult in c(0.25, 0.5, 1, 2)) {
  d <- e$vectors[, k] / D * sqrt(4 * tol / e$values[k]) * mult
  f0 <- f$obj$fn(p)
  cat(sprintf("step x %.2f: loss + %.4g  - %.4g (quadratic model %.4g)\n",
              mult, f$obj$fn(p + d) - f0, f$obj$fn(p - d) - f0,
              2 * tol * mult^2))
}
cat("direction loads:", paste0(nm, "=", round(e$vectors[, k], 3)), "\n")
ia <- match("a_(Intercept)", nm); ik <- match("k_(Intercept)", nm)
it <- match("theta_1", nm)
tg <- numeric(length(p)); tg[ia] <- p[ia]; tg[ik] <- -p[ik]; tg[it] <- 1
tg <- tg * D; tg <- tg / sqrt(sum(tg^2))
cat("|cos| with the ridge tangent:", abs(sum(tg * e$vectors[, k])), "\n")
