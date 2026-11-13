# Lane nlstop, items 4960 and 4902: the ridge-to-a-limit warning
# (nl_limit_message()). (1) the bernoulli Michaelis-Menten row (fresh set,
# seed 55) with xp perturbed by a few ulps, 20 draws: is the warning
# given on each rounding? (2) the objective at Vm, K times and divided by
# 1000. (3) absent cases: a body with a coefficient near 0 (b1 * exp(b2 *
# x), true b2 = 0), the gaussian Michaelis-Menten rows of the fresh set,
# and the sat body of the 450 rows (seed 101), each over 5 seeds.
# Usage: NLSTOP_ARM=after Rscript dev/nlstop-limit.R
source("dev/nlstop-common.R")
src <- parse("dev/nlstop-fresh.R")
keep <- vapply(src, function(e) {
  is.call(e) && identical(e[[1L]], as.name("<-")) &&
    as.character(e[[2L]])[1] %in% c("F", "pos", "B", "gen", "form")
}, NA)
for (e in src[keep]) eval(e)
lim <- function(r) any(attr(r, "wc") == "frmtmb_nl_limit")
d0 <- gen("bernoulli", "mm", 55 * 1000 + match("mm", names(B)))
fm <- form("bernoulli", "mm")
n_lim <- 0L
n_warn <- 0L
lls <- numeric(0)
for (k in 1:20) {
  d <- d0
  d$xp <- d$xp * (1 + 2e-15 * k)
  r <- tt(frm(fm, family = bernoulli(), data = d))
  n_lim <- n_lim + lim(r)
  n_warn <- n_warn + (length(attr(r, "w")) > 0L)
  lls <- c(lls, llv(r))
  cat(sprintf("ulp draw %d: ll %.6f code %s warn %s | %s\n", k, llv(r),
              r$opt$convergence, wsum(r),
              substr(paste(attr(r, "w"), collapse = " || "), 1, 120)))
}
cat(sprintf(paste("mm55 ulp draws: limit warning on %d of 20, any warning",
                  "on %d of 20, ll range %.3g\n"),
            n_lim, n_warn, diff(range(lls))))
r <- tt(frm(fm, family = bernoulli(), data = d0))
p <- r$opt$par
for (s in c(1e3, 1e-3)) {
  q <- p; q[1:2] <- q[1:2] * s
  cat(sprintf("Vm, K times %g: ll change %.3e\n", s,
              -(r$obj$fn(q) - r$opt$objective)))
}
# a coefficient near 0
for (seed in 1:5) {
  set.seed(seed)
  d <- data.frame(x = rnorm(300))
  d$y <- 2 + rnorm(300, 0, 0.5)
  r <- tt(frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), data = d))
  cat(sprintf("b2 near 0, seed %d: b2 %.3g limit %s warn %s\n", seed,
              r$estimates$beta[2], lim(r), wsum(r)))
}
for (fam in c("gaussian", "poisson")) for (seed in c(11, 22, 33, 44, 55)) {
  d <- gen(fam, "mm", seed * 1000 + match("mm", names(B)))
  r <- tt(frm(form(fam, "mm"), family = F[[fam]], data = d))
  cat(sprintf("%s mm seed %d: limit %s warn %s\n", fam, seed, lim(r),
              wsum(r)))
}
for (seed in 1:5) {
  set.seed(seed * 100 + 1)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  d$y <- 1.5 * exp(0.6 * d$x) + 0.3 * d$z + rnorm(n, 0, 0.5)
  r <- tt(frm(bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE),
              data = d))
  cat(sprintf("sat seed %d: limit %s warn %s\n", seed, lim(r), wsum(r)))
}
cat("DONE\n")
