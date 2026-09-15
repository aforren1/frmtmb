# What a random effect on gamma1 is anchored to, measured.
#
# eta = gamma0 + (gamma1 + u_c) x with x = log t. Rescale time by c and
# x' = x + log c, so the SAME data read in new units are
#   eta = (gamma0 - (gamma1 + u_c) log c) + (gamma1 + u_c) x'
# and the centre now carries a random INTERCEPT of -u_c log c beside
# its slope. A slope-only block is therefore a different model in
# different time units: it asserts that every centre has the same
# cumulative hazard at t = 1, and t = 1 is a unit rather than a fact.
#
# The log likelihood of a continuous density is not invariant to a
# rescale either: each EVENT row contributes log f and f' = f / c, so a
# refit in new units must be compared as ll + n_event log c. Censored
# rows contribute log S, which is invariant.
#
# SEED 20260910.
source("frailty-common.R")
suppressMessages({
  library(frmtmb); library(frmtmb.spline)
})

d <- g1_sim(20260910L)
nev <- sum(d$event)
bk <- range(log(d$time[d$event == 1L]))

fit_at <- function(mult, form) {
  dd <- d
  dd$time <- dd$time * mult
  fam <- royston_parmar(knots = numeric(0), bknots = bk + log(mult))
  suppressWarnings(frm(form, family = fam, data = dd, se = TRUE))
}

# Two CONTROLS ahead of the two arms, added in punch round 1.
#
# `no_block` has no centre term at all, so it is the same model in any
# time unit by construction and its spread MUST be flat. That is what
# licenses the `+ n_event log c` correction: without it a reader cannot
# separate a real anchoring from a Jacobian I got wrong, and the first
# draft of this measurement did get it wrong, by 2982 units.
#
# `fixed_slope` lets the slope vary by centre as FIXED effects with a
# shared intercept. If it moves too, the anchoring is a property of a
# varying slope with a common intercept rather than of the block being
# random, which is what the help page now says.
forms <- list(
  no_block = bf(time | cens(censored) ~ trt),
  fixed_slope = bf(time | cens(censored) ~ trt, gamma1 ~ centre),
  slope_only = bf(time | cens(censored) ~ trt, gamma1 ~ (1 | centre)),
  intercept_and_slope = bf(time | cens(censored) ~ trt + (1 | c | centre),
                           gamma1 ~ (1 | c | centre)))

mults <- c(1, 12, 365.25)
for (nm in names(forms)) {
  cat("\n==", nm, "==\n")
  adj <- c(); sds <- list()
  for (m in mults) {
    f <- fit_at(m, forms[[nm]])
    ll <- as.numeric(logLik(f)) + nev * log(m)
    adj <- c(adj, ll)
    vc <- frmtmb::VarCorr(f)
    sd <- unlist(lapply(vc, function(V) sqrt(diag(V))))
    if (!length(sd)) sd <- NA_real_
    sds[[as.character(m)]] <- sd
    if (m == 1 && length(vc)) V1 <- vc[[1L]]
    dg <- frmtmb::diagnose(f, quiet = TRUE)
    cat(sprintf(" mult %8.2f  ll+n_ev*log(mult) %.8f  sd %-28s  maxgrad %.2e\n",
                m, ll, paste(format(sd, digits = 6), collapse = " "),
                dg$max_grad))
  }
  cat(sprintf(" SPREAD over the three time units: %.4e log likelihood units\n",
              max(adj) - min(adj)))
  g1sd <- vapply(sds, function(s) s[[length(s)]], numeric(1))
  if (anyNA(g1sd)) {
    cat(" no variance component in this arm
")
    next
  }
  cat(sprintf(" sd(gamma1 | centre) over the three units: %s  relative spread %.3e\n",
              paste(format(g1sd, digits = 6), collapse = " "),
              (max(g1sd) - min(g1sd)) / mean(g1sd)))
}
cat("\nn_event =", nev, "\n")
