# Standard errors of a nested quadrature fit against the numerical
# Hessian of the plain-R reference integral (80 nodes per level).
#   Rscript dev/nestquad-se.R
# Data: the test's weighted Bernoulli design, seed 20261020.
.libPaths(c("C:/Users/adf44/source/r/wt-nestquad-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/nestquad-ref.R")
set.seed(20261020)
ng <- 12; nt <- 6; n <- ng * nt
d <- data.frame(ga = factor(rep(seq_len(ng), each = nt)),
                gb = factor(rep(c("i", "ii"), length.out = n)),
                x = rnorm(n), w = runif(n, 0.5, 2))
ab <- interaction(d$ga, d$gb, drop = TRUE)
d$eta <- 0.2 + 0.5 * d$x + rnorm(ng, 0, 0.7)[d$ga] +
  rnorm(nlevels(ab), 0, 1)[ab]
d$y <- rbinom(n, 1, plogis(d$eta))
fq <- frm(bf(y | weights(w) ~ 1 + x + (1 | ga/gb)), family = bernoulli(),
          data = d, quadrature = TRUE)
p <- fq$opt$par
print(p)
f <- function(q) {
  -ref_nested_gh("bernoulli", d$y, q[1] + q[2] * d$x, d$ga, ab, exp(q[3]),
                 exp(q[4]), w = d$w, k = 80)
}
h <- 1e-3
k <- length(p)
H <- matrix(0, k, k)
f0 <- f(p)
for (i in 1:k) for (j in i:k) {
  ei <- replace(numeric(k), i, h); ej <- replace(numeric(k), j, h)
  H[i, j] <- H[j, i] <- (f(p + ei + ej) - f(p + ei - ej) - f(p - ei + ej) +
                           f(p - ei - ej)) / (4 * h^2)
}
se_ref <- sqrt(diag(solve(H)))
sdr <- TMB::sdreport(fq$obj)
se_fit <- sqrt(diag(sdr$cov.fixed))
cat(sprintf("SE par %-6s fit %.6f ref %.6f ratio %.5f\n", names(p), se_fit,
            se_ref, se_fit / se_ref), sep = "")
cat(sprintf("OBJ fit %.10f ref %.10f\n", fq$opt$objective, f0))
