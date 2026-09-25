# The conditional draw of dropped columns of a partially kept term,
# against the exact conditional Gaussian, internal and end to end.
#   Rscript dev/simnewdata-review/rv-conditional.R > .../log/conditional.txt
source("dev/simnewdata-review/rv-prelude.R")
set.seed(31)
n_g <- 30
d <- data.frame(x1 = rnorm(n_g * 25), x2 = rnorm(n_g * 25),
                g = factor(rep(seq_len(n_g), 25)))
Vt <- matrix(c(1, .6, -.5, .6, .8, -.2, -.5, -.2, .6), 3)
U <- matrix(rnorm(n_g * 3), n_g) %*% chol(Vt)
d$y <- 1 + d$x1 - d$x2 + U[d$g, 1] + U[d$g, 2] * d$x1 + U[d$g, 3] * d$x2 +
  rnorm(nrow(d), 0, 0.5)
fit <- frm(bf(y ~ x1 + x2 + (1 + x1 + x2 | g)), data = d)
V <- varcorr_matrices(fit)[[1]]
cat("V (fitted):\n"); print(round(V, 4))
bk <- fit$frame$re_blocks[[1]]
B0 <- matrix(fit$estimates$b[bk$b_idx], nrow = 3)
# cross-check the b layout against ranef()
re <- ranef(fit)
cat("layout check, level 1 b vs ranef:", format(B0[, 1], digits = 5), "|",
    format(unlist(re[[1]][1, ]), digits = 5), "\n")

cond_exact <- function(V, K, bK) {
  D <- setdiff(1:3, K)
  A <- V[D, K, drop = FALSE] %*% solve(V[K, K, drop = FALSE])
  list(D = D, m = as.vector(A %*% bK),
       C = V[D, D, drop = FALSE] - A %*% V[K, D, drop = FALSE])
}

R <- 20000
internal <- function(rf, K, lev = 4) {
  plan <- frmtmb:::sim_re_plan(fit, rf)
  set.seed(1)
  dr <- t(replicate(R, {
    matrix(frmtmb:::sim_draw_b(fit, plan)[bk$b_idx], nrow = 3)[, lev]
  }))
  ex <- cond_exact(V, K, B0[K, lev])
  kept_ok <- all(dr[, K] == rep(B0[K, lev], each = R))
  mz <- (colMeans(dr[, ex$D, drop = FALSE]) - ex$m) /
    sqrt(diag(ex$C) / R)
  E <- stats::cov(dr[, ex$D, drop = FALSE])
  sdc <- sqrt((outer(diag(ex$C), diag(ex$C)) + ex$C^2) / R)
  cat(sprintf("internal %-18s kept cols %s fixed=%s  mean max|z| %.2f  cov max|z| %.2f\n",
              deparse1(rf), paste(K, collapse = ","), kept_ok,
              max(abs(mz)), max(abs((E - ex$C) / sdc))))
}
internal(~ (1 | g), 1)
internal(~ (1 + x1 | g), 1:2)
internal(~ (0 + x2 | g), 3)
internal(~ (x2 + 1 | g), c(1, 3))

# end to end: one fitted level, three rows
nd <- data.frame(x1 = c(0, 1, -1), x2 = c(0, -1, 2), g = factor("4"))
Z <- cbind(1, nd$x1, nd$x2)
beta <- fixef(fit); if (is.matrix(beta)) beta <- beta[, "Estimate"]
e2e <- function(rf, K) {
  s <- as.matrix(simulate(fit, nsim = R, seed = 7, re_formula = rf,
                          newdata = nd))
  ex <- cond_exact(V, K, B0[K, 4])
  bmean <- numeric(3); bmean[K] <- B0[K, 4]; bmean[ex$D] <- ex$m
  mu <- as.vector(cbind(1, nd$x1, nd$x2) %*% beta + Z %*% bmean)
  S <- Z[, ex$D, drop = FALSE] %*% ex$C %*% t(Z[, ex$D, drop = FALSE]) +
    diag(sigma(fit)^2, 3)
  E <- stats::cov(t(s))
  sdc <- sqrt((outer(diag(S), diag(S)) + S^2) / R)
  cat(sprintf("e2e      %-18s mean max|z| %.2f  cov max|z| %.2f\n",
              deparse1(rf), max(abs((rowMeans(s) - mu) / sqrt(diag(S) / R))),
              max(abs((E - S) / sdc))))
  # the marginal alternative, to show the probe can tell them apart
  Sm <- Z[, ex$D, drop = FALSE] %*% V[ex$D, ex$D] %*%
    t(Z[, ex$D, drop = FALSE]) + diag(sigma(fit)^2, 3)
  cat(sprintf("         (marginal law would give cov max|z| %.2f, mean max|z| %.2f)\n",
              max(abs((E - Sm) / sdc)),
              max(abs((rowMeans(s) - as.vector(cbind(1, nd$x1, nd$x2) %*% beta +
                 Z[, K, drop = FALSE] %*% B0[K, 4])) / sqrt(diag(S) / R)))))
}
e2e(~ (1 | g), 1)
e2e(~ (1 + x1 | g), 1:2)

# the same with a second, uncorrelated term, so draw_b() also runs
d$h <- factor(rep(1:6, length.out = nrow(d)))
d$y <- d$y + rnorm(6, 0, 0.7)[d$h]
fit2 <- frm(bf(y ~ x1 + x2 + (1 + x1 + x2 | g) + (1 | h)), data = d)
p2 <- frmtmb:::sim_re_plan(fit2, ~ (1 | g))
cat("two-term plan: redraw", length(p2$redraw), "cond", length(p2$cond),
    "\n")

# refusals: every structure without a conditional law
refuse <- function(label, f, dat, rf) {
  ft <- tryCatch(suppressWarnings(frm(f, data = dat)),
                 error = function(e) e)
  if (inherits(ft, "error")) {
    cat(sprintf("%-10s FIT ERROR: %s\n", label, conditionMessage(ft)))
    return(invisible())
  }
  cs <- vapply(ft$frame$re_blocks, `[[`, "", "covstruct")
  cat(sprintf("%-10s covstructs %s: %s\n", label, paste(cs, collapse = ","),
              substr(try_msg(simulate(ft, nsim = 2, seed = 1,
                                      re_formula = rf)), 1, 160)))
}
A <- diag(n_g); A[A == 0] <- 0.2
dimnames(A) <- list(levels(d$g), levels(d$g))
refuse("gr_cov", bf(y ~ x1 + (1 + x1 | gr(g, cov = A))), d, ~ (1 | g))
refuse("student", bf(y ~ x1 + (1 + x1 | gr(g, dist = "student"))), d,
       ~ (1 | g))
refuse("rr", bf(y ~ x1 + rr(1 + x1 + x2 | g, d = 1)), d, ~ (1 | g))
refuse("cs", bf(y ~ x1 + cs(1 + x1 | g)), d, ~ (1 | g))
refuse("diag", bf(y ~ x1 + diag(1 + x1 | g)), d, ~ (1 | g))
refuse("ar1", bf(y ~ x1 + ar1(0 + h | g)), d, ~ (0 + h1 | g))
