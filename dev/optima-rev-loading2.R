# Reviewer of lane optima, claim 7: the lost test's property rebuilt
# without mo(). A gaussian fit's own objective, with a Hessian built to
# have one removed direction that is mostly `z` and loads 0.022 on `x`
# (the seed-12 mo() figure). se_tier3() must remove z and keep x, and
# the null basis it hands predictions must have a zero in x's row.
# Seen to fail: the same call on se_tier3() without the line that
# zeroes a kept parameter's loading.
#   Rscript dev/optima-rev-loading2.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(x = rnorm(200), z = rnorm(200))
d$y <- 1 + 0.5 * d$x + 0.3 * d$z + rnorm(200)
fit <- frm(y ~ x + z, data = d)
nm <- frmtmb:::outer_par_names(fit)
p <- fit$opt$par
Hr <- fit$obj$he(p)
D <- sqrt(diag(Hr))
np <- length(p)
ix <- match("x", sub("^.*_", "", nm))
iz <- match("z", sub("^.*_", "", nm))
cat("params:", nm, "| x at", ix, "z at", iz, "\n")
v1 <- numeric(np)
v1[iz] <- 1
v1[ix] <- 0.022
v1 <- v1 / sqrt(sum(v1^2))
Q <- qr.Q(qr(cbind(v1, diag(np))))[, seq_len(np)]
if (sum(Q[, 1] * v1) < 0) Q[, 1] <- -Q[, 1]
for (lam in c(-0.5, 0)) {
  S <- Q %*% diag(c(lam, rep(1, np - 1))) %*% t(Q)
  H <- S * outer(D, D)
  E <- matrix(0, np, np)
  t3 <- frmtmb:::se_tier3
  mut <- t3
  b <- deparse(body(mut))
  hit <- grep("Vr[mark == \"\", ] <- 0", b, fixed = TRUE)
  body(mut) <- parse(text = b[-hit])[[1]]
  environment(mut) <- asNamespace("frmtmb")
  a <- t3(fit, H, E, p, exact = TRUE)
  m <- mut(fit, H, E, p, exact = TRUE)
  cat(sprintf(paste0("lambda %4.1f: lost %s | null row x: lane %s, ",
                     "mutant %s\n"), lam,
              paste(names(a$lost), a$lost, sep = "=", collapse = ","),
              paste(format(a$null[ix, ], digits = 3), collapse = ","),
              paste(format(m$null[ix, ], digits = 3), collapse = ",")))
}
