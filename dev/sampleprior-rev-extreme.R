# Reviewer of lane sampleprior, fix round 1: the extreme-draw check in
# draws_cov_values() against planted region defects in the fast path.
# Model and seeds as dev/sampleprior-rev-names.R.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
ns <- asNamespace("frmtmb.sample")
set.seed(20261118)
G <- 12
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(rep(seq_len(G), each = n / G)),
                h = factor(rep(1:4, length.out = n)))
d$g2 <- factor(rep(1:3, length.out = n))
S <- matrix(c(1, .5, -.3, .5, .8, .2, -.3, .2, .6), 3)
B <- MASS::mvrnorm(G, c(0, 0, 0), S)
d$y <- 1 + B[d$g, 1] + (1 + B[d$g, 2]) * d$x + (0.5 + B[d$g, 3]) * d$z +
  rnorm(4)[d$h] + rnorm(G * 3)[as.integer(interaction(d$g, d$g2))] +
  rnorm(n, 0, 0.7)
fit <- frm(y ~ x + z + (1 + x + z | g) + (1 | h) + (1 | g:g2), data = d)
ds <- suppressMessages(suppressWarnings(frm_sample(
  fit, chains = 2, iter = 600, refresh = 0, seed = 21)))
m <- ns$draws_internal_matrix(ds)
dv <- ns$draws_natural_cols(ds$fit)$derived
slow <- dv; slow$fast <- NULL
ref <- ns$draws_cov_values(m, ds$fit, slow)
nm <- dv$fast$names
plant <- function(lab, bad) {
  dm <- dv
  f0 <- dv$fast$vals_at
  dm$fast$vals_at <- function(theta) {
    v <- f0(theta)
    names(v) <- nm
    bad(v)
  }
  got <- ns$draws_cov_values(m, ds$fit, dm)
  err <- max(abs(got - ref))
  cat(sprintf("%-52s rows hit %4d | max |got - slow| %.3g -> %s\n", lab,
              sum(rowSums(abs(got - ref) > 0) > 0), err,
              if (err == 0) "caught (slow path)" else "MISSED"))
}
sdx <- ref[, "sd_g__x"]
cxz <- ref[, "cor_g__x__z"]
cat("sd_g__x range", range(sdx), " cor_g__x__z range", range(cxz), "\n")
plant("sd_g__x off by 1e-6 where sd > its 90th pct",
      function(v) { if (v["sd_g__x"] > quantile(sdx, .9)) v["sd_g__x"] <- v["sd_g__x"] + 1e-6; unname(v) })
plant("cor_g__x__z off by 1e-6 where cor > its 99th pct",
      function(v) { if (v["cor_g__x__z"] > quantile(cxz, .99)) v["cor_g__x__z"] <- v["cor_g__x__z"] + 1e-6; unname(v) })
q <- quantile(cxz, c(.45, .55))
plant("cor_g__x__z off by 1e-6 inside its 45-55th pct",
      function(v) { if (v["cor_g__x__z"] > q[1] && v["cor_g__x__z"] < q[2]) v["cor_g__x__z"] <- v["cor_g__x__z"] + 1e-6; unname(v) })
