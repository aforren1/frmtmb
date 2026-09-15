# Does rp_floored() SEE the centre deviations?
#
# The gamma1 recovery row reports 0 non-monotone rows over 60 fits. A
# check that only ever reads the population gamma1 would report the
# same 0, so the claim needs the case where a centre's own slope is the
# thing that goes non-positive.
#
# The construction: a small population shape with a wide centre spread
# and few subjects per centre, so the conditional modes are noisy
# enough to cross zero. The truth is reject sampled to keep every true
# shape positive, since t^(negative) is not a survival distribution;
# what is being tested is the CHECK, not the estimator.
#
# SEEDS 20260910 upward.
source("frailty-common.R")
suppressMessages({
  library(frmtmb); library(frmtmb.spline)
})

sim_thin <- function(seed, n = 400L, n_centre = 40L, sd_u = 0.35,
                     gamma1 = 0.6, scale = 5, beta = 0.6,
                     p_cens = 0.4, floor_true = 0.05) {
  set.seed(seed)
  centre <- rep(seq_len(n_centre), length.out = n)
  repeat {
    u <- stats::rnorm(n_centre, 0, sd_u)
    if (all(gamma1 + u > floor_true)) break
  }
  trt <- stats::rbinom(n, 1L, 0.5)
  g0 <- -gamma1 * log(scale)
  tt <- exp((log(-log(stats::runif(n))) - g0 - beta * trt) /
              (gamma1 + u[centre]))
  tau <- unname(stats::quantile(tt, 1 - p_cens))
  ev <- as.integer(tt <= tau)
  d <- data.frame(time = pmin(tt, tau), event = ev, censored = 1L - ev,
                  trt = trt, centre = factor(centre))
  attr(d, "truth") <- list(shape = gamma1 + u)
  d
}

hits <- 0L
for (s in 20260910L + 0:14) {
  d <- sim_thin(s)
  bk <- range(log(d$time[d$event == 1L]))
  f <- try(suppressWarnings(frm(
    bf(time | cens(censored) ~ trt, gamma1 ~ (1 | centre)),
    family = royston_parmar(knots = numeric(0), bknots = bk),
    data = d, se = TRUE)), silent = TRUE)
  if (inherits(f, "try-error")) {
    cat(s, "fit failed\n")
    next
  }
  g1 <- unname(fixef(f)$gamma1[["(Intercept)"]])
  sh <- g1 + as.numeric(frmtmb::ranef(f)[["centre"]])
  rep_ <- rp_floored(f, action = "report")
  err <- try(rp_floored(f), silent = TRUE)
  cat(sprintf(
    "seed %d  pop gamma1 %7.4f  min centre slope %8.4f  n_nonmono %3d  refuses %s\n",
    s, g1, min(sh), rep_$n_nonmonotone,
    inherits(err, "try-error")))
  if (rep_$n_nonmonotone > 0L) {
    hits <- hits + 1L
    # the rows named must belong to the centres whose slope is negative
    rows <- attr(rep_, "rows")[["nonmonotone"]]
    cc <- unique(as.integer(as.character(d$centre[rows])))
    neg <- which(sh <= 0)
    cat("   rows come from centres ", paste(sort(cc), collapse = " "),
        "; centres with a non-positive fitted slope are ",
        paste(neg, collapse = " "), "\n", sep = "")
    cat("   population gamma1 alone is ",
        if (g1 > 0) "POSITIVE, so a population-only check reports 0"
        else "non-positive too, which would not separate the two",
        "\n", sep = "")
  }
}
cat("\nseeds with the floor firing:", hits, "of 15\n")
