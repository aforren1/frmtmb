# Reviewer 2, items 3/4: left- and interval-censored rows under the
# contaminant, against a hand-written mixture with RWiener. No fitting.
source("dev/phase3b-review2/r2-prelude.R"); r2_lib("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam); library(RWiener)})
set.seed(31)
d <- ddm_simulate(250, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
hit <- runif(250) < 0.1
d$rt[hit] <- runif(sum(hit), 0.2, 3); d$upper[hit] <- rbinom(sum(hit), 1, 0.5)
d$code <- 0L; d$y2 <- d$rt
l <- d$rt < 0.5; d$code[l] <- -1L; d$rt[l] <- 0.5
pk <- which(d$code == 0L)[seq(2, sum(d$code == 0L), by = 3)]
lo <- floor(d$rt[pk] * 10) / 10
d$code[pk] <- 2L; d$y2[pk] <- lo + 0.1; d$rt[pk] <- pmax(lo, 0.5)
r <- d$code == 0L & d$rt > 2; d$code[r] <- 1L; d$rt[r] <- 2
print(table(d$code))
cr <- c(0.1, 3)
o <- frm(bf(rt | dec(upper) + cens(code, y2) ~ 1, bs ~ 1, ndt = 0.3, lambda ~ 1, bias = 0.5),
         family = wiener(contaminant = TRUE, contaminant_range = cr, max_ndt = 0.45),
         data = d, dry_run = "objective")
cat("parameters:", names(o$obj$par), "\n")
G <- function(y) pmin(pmax((y - cr[1]) / diff(cr), 0), 1)
hand <- function(p) {
  v <- p[1]; a <- exp(p[2]); lam <- plogis(p[3]); t0 <- 0.3; w <- 0.5
  resp <- ifelse(d$upper == 1, "upper", "lower")
  ll <- vapply(seq_len(nrow(d)), function(i) {
    y <- d$rt[i]
    switch(as.character(d$code[i]),
      "0" = log((1 - lam) * dwiener(y, a, t0, w, v, resp = resp[i]) + lam / (2 * diff(cr))),
      "1" = log((1 - lam) * (1 - pwiener(y, a, t0, w, v, resp = "both")) + lam * (1 - G(y))),
      "-1" = log((1 - lam) * pwiener(y, a, t0, w, v, resp = resp[i]) + lam * G(y) / 2),
      "2" = log((1 - lam) * (pwiener(d$y2[i], a, t0, w, v, resp = resp[i]) -
                              pwiener(y, a, t0, w, v, resp = resp[i])) +
                  lam * (G(d$y2[i]) - G(y)) / 2))
  }, 0)
  -sum(ll)
}
for (p in list(o$obj$par, c(0.5, 0.2, -2), c(1.2, 0.5, -4))) {
  cat(sprintf("lane %.10f  hand %.10f  rel %.2e\n", o$obj$fn(p), hand(p),
              abs(o$obj$fn(p) - hand(p)) / abs(hand(p))))
}
