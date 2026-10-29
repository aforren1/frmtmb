## Review rescorpp: is predict_rescor_draw() MVN(mu, D R D) / MVT with
## row-varying sigma? Direct calls on ML fits, then draw/seed semantics
## on draws. Rscript dev/review-rescorpp-math.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
library(MASS)

set.seed(101)
n <- 150
x <- rnorm(n)
R0 <- matrix(c(1, .8, -.5, .8, 1, -.3, -.5, -.3, 1), 3)
E <- mvrnorm(n, rep(0, 3), R0)
dd <- data.frame(x = x, y1 = 1 + x + exp(1.5 * x) * E[, 1],
                 y2 = 100 * E[, 2], y3 = 0.001 * E[, 3])
fg <- quiet(frm(mvbf(bf(y1 ~ x, sigma ~ x), bf(y2 ~ x), bf(y3 ~ x)) +
                  set_rescor(TRUE), family = gaussian(), data = dd))
R <- rescor_matrix(fg)
cat("rescor_matrix:\n"); print(round(R, 4))
dp <- eval_dpars(fg)
s1 <- dp$y1$sigma; s2 <- dp$y2$sigma; s3 <- dp$y3$sigma
cat("sigma ranges: y1", signif(range(s1), 3), " y2", signif(range(s2), 3),
    " y3", signif(range(s3), 3), "\n")
cat("is predict_rescor_draw exported:",
    "predict_rescor_draw" %in% getNamespaceExports("frmtmb"), "\n")

## M replicates of every row
M <- 20000
ok <- rep(TRUE, n)
set.seed(7)
Y <- array(NA_real_, c(M, n, 3))
for (m in seq_len(M)) {
  ys <- predict_rescor_draw(fg, fg$spec$responses, dp, ok)
  Y[m, , 1] <- ys$y1; Y[m, , 2] <- ys$y2; Y[m, , 3] <- ys$y3
}
## per row: whiten by chol(D_i R D_i); the result must be iid N(0, I)
mu <- cbind(dp$y1$mu, dp$y2$mu, dp$y3$mu)
sg <- cbind(s1, rep_len(s2, n), rep_len(s3, n))
ks_p <- numeric(0); offcor <- numeric(0); vz <- numeric(0)
pick <- unique(c(which.min(s1), which.max(s1), round(seq(1, n, length.out = 8))))
for (i in pick) {
  Si <- diag(sg[i, ]) %*% R %*% diag(sg[i, ])
  Ui <- sweep(Y[, i, ], 2, mu[i, ]) %*% solve(chol(Si))
  for (k in 1:3) {
    ks_p <- c(ks_p, suppressWarnings(ks.test(Ui[, k], "pnorm")$p.value))
    vz <- c(vz, (var(Ui[, k]) - 1) / sqrt(2 / (M - 1)))
  }
  cc <- cor(Ui)
  offcor <- c(offcor, cc[upper.tri(cc)] * sqrt(M))
  # Mahalanobis distance against chi^2_3
  q <- rowSums(Ui^2)
  ks_p <- c(ks_p, suppressWarnings(ks.test(q, "pchisq", 3)$p.value))
}
cat("rows checked:", length(pick), " sigma_y1 at them:",
    signif(s1[pick], 3), "\n")
cat("gaussian whitened: KS p min", signif(min(ks_p), 3), "of", length(ks_p),
    "tests; var z range", round(range(vz), 2),
    "; off-diag cor*sqrt(M) range", round(range(offcor), 2), "\n")
## and the negative control: drop the Cholesky factor (independent draws)
Zi <- matrix(rnorm(M * 3), M) # independent N(0,1) per response
i <- pick[1]
Si <- diag(sg[i, ]) %*% R %*% diag(sg[i, ])
Ui0 <- (Zi * rep(sg[i, ], each = M)) %*% solve(chol(Si))
cat("control (independent) row", i, ": KS p of Mahalanobis",
    signif(suppressWarnings(ks.test(rowSums(Ui0^2), "pchisq", 3)$p.value), 3),
    "\n")

## Student-t: shared W per row
set.seed(102)
Et <- mvrnorm(n, rep(0, 3), R0) / sqrt(rchisq(n, 5) / 5)
dt3 <- data.frame(x = x, y1 = 1 + x + Et[, 1], y2 = 3 * Et[, 2],
                  y3 = 0.2 * Et[, 3])
ft <- quiet(frm(mvbf(bf(y1 ~ x, sigma ~ x), bf(y2 ~ x), bf(y3 ~ x)) +
                  set_rescor(TRUE), family = student(), data = dt3))
Rt <- rescor_matrix(ft)
dpt <- eval_dpars(ft)
nu <- frmtmb:::rescor_shared_nu(ft)
cat("\nstudent nu:", nu, " nu in dpars y1:", format(dpt$y1$nu)[1],
    " sigma_y1 range", signif(range(dpt$y1$sigma), 3), "\n")
set.seed(8)
Yt <- array(NA_real_, c(M, n, 3))
for (m in seq_len(M)) {
  ys <- predict_rescor_draw(ft, ft$spec$responses, dpt, ok)
  Yt[m, , 1] <- ys$y1; Yt[m, , 2] <- ys$y2; Yt[m, , 3] <- ys$y3
}
mut <- cbind(dpt$y1$mu, dpt$y2$mu, dpt$y3$mu)
sgt <- cbind(rep_len(dpt$y1$sigma, n), rep_len(dpt$y2$sigma, n),
             rep_len(dpt$y3$sigma, n))
ksm <- numeric(0); ksF <- numeric(0)
for (i in pick) {
  for (k in 1:3) {
    zk <- (Yt[, i, k] - mut[i, k]) / sgt[i, k]
    ksm <- c(ksm, suppressWarnings(ks.test(zk, "pt", nu)$p.value))
  }
  Si <- diag(sgt[i, ]) %*% Rt %*% diag(sgt[i, ])
  Ui <- sweep(Yt[, i, ], 2, mut[i, ]) %*% solve(chol(Si))
  ksF <- c(ksF, suppressWarnings(ks.test(rowSums(Ui^2) / 3, "pf", 3, nu)$p.value))
}
cat("student margins vs pt(nu): KS p min", signif(min(ksm), 3), "of",
    length(ksm), "\n")
cat("student Mahalanobis/3 vs F(3, nu): KS p min", signif(min(ksF), 3),
    "of", length(ksF), "\n")
## against brms's own rmulti_student_t at one row (two-sample KS)
i <- pick[2]
Si <- diag(sgt[i, ]) %*% Rt %*% diag(sgt[i, ])
set.seed(9)
Bb <- brms::rmulti_student_t(M, df = nu, mu = mut[i, ], Sigma = Si)
qa <- mahalanobis(Yt[, i, ], mut[i, ], Si)
qb <- mahalanobis(Bb, mut[i, ], Si)
cat("row", i, "two-sample KS (frmtmb vs brms rmulti_student_t) on",
    "Mahalanobis p", signif(suppressWarnings(ks.test(qa, qb)$p.value), 3),
    "; on y1 margin p", signif(suppressWarnings(
      ks.test(Yt[, i, 1], Bb[, 1])$p.value), 3), "\n")
## negative control: separate W per response breaks the F law
Wsep <- sqrt(matrix(rchisq(M * 3, nu), M) / nu)
Zs <- (matrix(rnorm(M * 3), M) %*% chol(Rt)) / Wsep
cat("control (separate W): KS p of Mahalanobis/3 vs F",
    signif(suppressWarnings(ks.test(rowSums((Zs %*% solve(chol(Rt)))^2) / 3,
                                    "pf", 3, nu)$p.value), 3), "\n")

## subset of responses through predict_rescor_draw: c("y3", "y1")
rs <- ft$spec$responses[c("y3", "y1")]
set.seed(10)
Y31 <- t(vapply(seq_len(5000), function(m) {
  ys <- predict_rescor_draw(ft, rs, dpt[c("y3", "y1")], ok)
  c(ys$y3[1], ys$y1[1])
}, numeric(2)))
cat("subset (y3, y1) row 1 Kendall tau", round(cor(Y31[, 1], Y31[, 2],
    method = "kendall"), 4), " target 2/pi asin(r)",
    round(2 / pi * asin(Rt["y1", "y3"]), 4), "\n")
cat("names returned:", names(predict_rescor_draw(ft, rs, dpt[c("y3", "y1")], ok)), "\n")

## ok partially FALSE: rows dropped, others still correct
ok2 <- rep(c(TRUE, FALSE), length.out = n)
ys <- predict_rescor_draw(fg, fg$spec$responses, dp, ok2)
cat("ok half: lengths", lengths(ys), "\n")
