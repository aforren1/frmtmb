# Reviewer of lane optima, task 1: the prior a flat (Lebesgue) density
# on the free coordinates implies on the mo() simplex, softmax (base)
# against the stereographic chart (lane).
#   Rscript dev/optima-rev-sample1-implied.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
ms_lane <- utils::getFromNamespace("mo_simplex", "frmtmb")
basis <- utils::getFromNamespace("mo_sphere_basis", "frmtmb")
cat("frmtmb from", find.package("frmtmb"), "\n")

# vectorized copies; checked against the package's mo_simplex below
chart <- function(Z, D) {
  B <- basis(D); P <- -1 / sqrt(D)
  r2 <- rowSums(Z * Z)
  U <- (2 * Z %*% t(B) + (r2 - 1) * P) / (r2 + 1)
  W <- U * U
  W / rowSums(W)
}
softmax <- function(Z) {
  E <- exp(cbind(0, Z))
  E / rowSums(E)
}
set.seed(11)
Zt <- matrix(rnorm(30, sd = 3), 10)
stopifnot(max(abs(chart(Zt, 4) - t(apply(Zt, 1, ms_lane)))) < 1e-12)
cat("vectorized chart matches frmtmb:::mo_simplex (D = 4, 10 points)\n")

# per-sheet preimage of w under the chart: u = s * sqrt(w)
preimage <- function(w, s, D) {
  u <- s * sqrt(w)
  as.numeric(crossprod(basis(D), u)) / (1 + sum(u) / sqrt(D))
}
# reviewer's density, with constants: 2^(1-D) prod(w)^(-1/2)
#   sum_s 1{z_s in box} (1 + sum(s sqrt w)/sqrt D)^(-(D-1)),
# w.r.t. dw_1..dw_{D-1}; the total mass equals the box volume (2R)^(D-1)
signs <- function(D) as.matrix(expand.grid(rep(list(c(-1, 1)), D)))
dens_chart <- function(w, D, R = Inf) {
  S <- signs(D); tot <- 0
  for (k in seq_len(nrow(S))) {
    s <- S[k, ]
    z <- preimage(w, s, D)
    if (isTRUE(all(abs(z) <= R))) {
      tot <- tot + (1 + sum(s * sqrt(w)) / sqrt(D))^(-(D - 1))
    }
  }
  2^(1 - D) * prod(w)^(-1 / 2) * tot
}
# softmax: zeta_k = log(w_{k+1} / w_1); |Jacobian| = 1 / prod(w)
dens_soft <- function(w, D, R = Inf) {
  z <- log(w[-1] / w[1])
  if (isTRUE(all(abs(z) <= R))) 1 / prod(w) else 0
}

# ---- verification, D = 3, 2D binning on (w1, w2) ----
verify <- function(map, dens, R, N = 2e6, nb = 10, D = 3, sub = 40) {
  Z <- matrix(runif(N * (D - 1), -R, R), N)
  W <- map(Z)
  vol <- (2 * R)^(D - 1)
  h <- 1 / nb
  i1 <- pmin(floor(W[, 1] / h), nb - 1); i2 <- pmin(floor(W[, 2] / h), nb - 1)
  emp <- table(factor(i1, 0:(nb - 1)), factor(i2, 0:(nb - 1)))
  out <- NULL
  for (a in 0:(nb - 1)) for (b in 0:(nb - 1)) {
    if (a + b > nb - 2) next          # keep bins wholly inside the simplex
    if (a == 0 || b == 0 || a + b == nb - 2) next  # skip face bins
    g <- (seq_len(sub) - 0.5) / sub * h
    pts <- expand.grid(w1 = a * h + g, w2 = b * h + g)
    pd <- mean(apply(pts, 1, function(p) {
      dens(c(p[1], p[2], 1 - p[1] - p[2]), D, R)
    })) * h * h / vol              # predicted probability of the bin
    ob <- emp[a + 1, b + 1] / N
    se <- sqrt(pd * (1 - pd) / N)
    out <- rbind(out, data.frame(a = a, b = b, pred = pd, obs = ob,
                                 z = (ob - pd) / se))
  }
  out
}
set.seed(2026)
for (R in c(5, 20)) {
  vc <- verify(function(Z) chart(Z, 3), dens_chart, R)
  vs <- verify(softmax, dens_soft, R)
  cat(sprintf(paste0("VERIFY R = %g: chart  %d interior bins, max |obs-pred|",
                     "/pred = %.4f, max |z| = %.2f, sum z^2 = %.1f\n"),
              R, nrow(vc), max(abs(vc$obs - vc$pred) / vc$pred),
              max(abs(vc$z)), sum(vc$z^2)))
  cat(sprintf(paste0("VERIFY R = %g: softmax %d interior bins, max |obs-pred|",
                     "/pred = %.4f, max |z| = %.2f, sum z^2 = %.1f\n"),
              R, nrow(vs), max(abs(vs$obs - vs$pred) / vs$pred),
              max(abs(vs$z)), sum(vs$z^2)))
  cat("  chart worst bins (a, b index w1, w2 in steps of 0.1):\n")
  print(head(vc[order(-abs(vc$z)), ], 4), digits = 4)
}
# a deliberately wrong alternative to show the test has power: the chart
# density with only the all-positive sheet
dens_wrong <- function(w, D, R = Inf) {
  2^(1 - D) * prod(w)^(-1 / 2) * 2^D *
    (1 + sum(sqrt(w)) / sqrt(D))^(-(D - 1))
}
set.seed(2026)
vw <- verify(function(Z) chart(Z, 3), dens_wrong, 5)
cat(sprintf("POWER: one-sheet wrong formula, R = 5: max |z| = %.1f\n",
            max(abs(vw$z))))

# D = 2 closed form check: w1 density at a few points vs histogram
set.seed(7)
for (R in c(5, 20)) {
  z <- matrix(runif(4e6, -R, R))
  w1 <- chart(z, 2)[, 1]
  br <- seq(0, 1, by = 0.05)
  hc <- hist(w1, br, plot = FALSE)$counts / length(w1) / 0.05 * 2 * R
  mid <- br[-1] - 0.025
  pc <- sapply(mid, function(m) dens_chart(c(m, 1 - m), 2, R))
  cat(sprintf("D = 2 R = %g: w1 density obs/pred at mids %s\n", R,
              paste(sprintf("%.3f", (hc / pc)[c(2, 5, 8, 10, 11, 14, 19)]),
                    collapse = " ")))
}

# ---- where the infinite mass sits: integrals over shrinking sets ----
cat("\nD = 3 density at chosen points (no box):\n")
for (w in list(c(1, 1, 1) / 3, c(.34, .33, .33), c(.4, .3, .3),
               c(.6, .3, .1), c(.98, .01, .01), c(.5, .49, .01))) {
  cat(sprintf("  w = (%s): chart %.4g  softmax %.4g\n",
              paste(format(w, digits = 2), collapse = ", "),
              dens_chart(w, 3), dens_soft(w, 3)))
}

# ---- box-uniform fractions ----
cat("\nFractions of box-uniform zeta draws (N = 1e6 each):\n")
cat("  near_bary: ||w - 1/D||_2 < 0.05;  near_face: min(w) < 0.01\n")
set.seed(99)
for (D in c(2, 3, 4)) for (R in c(5, 20, 100)) {
  Z <- matrix(runif(1e6 * (D - 1), -R, R), ncol = D - 1)
  for (nm in c("softmax", "chart")) {
    W <- if (nm == "softmax") softmax(Z) else chart(Z, D)
    db <- sqrt(rowSums((W - 1 / D)^2))
    cat(sprintf("  D = %d R = %3g %-7s near_bary %.4f near_face %.4f  median min(w) %.3g\n",
                D, R, nm, mean(db < 0.05), mean(do.call(pmin, as.data.frame(W)) < 0.01),
                median(do.call(pmin, as.data.frame(W)))))
  }
}
