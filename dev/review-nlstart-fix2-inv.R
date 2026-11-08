# Review of lane nlstart, second short check: scale invariance. The same
# data as review-nlstart-fix2-small.R; the body on x and on xs = x * s
# (s = 1e-4, 1e-6) is one model with b2 rescaled, so a default fit that
# does not depend on the covariate's scale reads the same logLik.
# Usage: REV_ARM=after|before Rscript dev/review-nlstart-fix2-inv.R
source("dev/review-nlstart-common.R")
F <- list(gaussian = gaussian(), poisson = poisson(), bernoulli = bernoulli(),
          cumulative = cumulative())
cat("ROWHEAD body fam seed ll_x ll_1e-4 ll_1e-6 w_x w_4 w_6\n")
for (bn in c("exp", "sat")) for (fam in names(F)) for (seed in 1:5) {
  set.seed(seed * 37 + 11)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), xp = round(runif(n, 0, 2), 3))
  lat <- 1.5 * exp(0.6 * d$x)
  d$y <- switch(fam, gaussian = lat + rnorm(n, 0, 0.3),
                poisson = rpois(n, pmax(lat, 0.05)),
                bernoulli = rbinom(n, 1, plogis(lat - 1)),
                cumulative = ordered(cut(lat + rlogis(n),
                                         c(-Inf, 0.3, 1, 2, Inf), labels = 1:4)))
  v <- c(); w <- c()
  for (s in c(1, 1e-4, 1e-6)) {
    d$xs <- d$x * s
    f <- if (bn == "exp") "b1 * exp(b2 * xs)" else "b1 * (1 - exp(-b2 * xs))"
    r <- tt(frm(bf(as.formula(paste("y ~", f)), b1 + b2 ~ 1, nl = TRUE),
                family = F[[fam]], data = d))
    v <- c(v, llv(r)); w <- c(w, length(attr(r, "w")))
  }
  cat(sprintf("ROW %s %s %d %.6f %.6f %.6f %d %d %d\n", bn, fam, seed, v[1],
              v[2], v[3], w[1], w[2], w[3]))
  flush.console()
}
cat("DONE\n")
