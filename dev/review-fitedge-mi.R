# Review of lane fitedge, item 4: where latres's mi() refusal fires. Args: arm.
# Data: set.seed(73000 + seed), 10 groups of 6 times, an AR(1) latent.
a <- commandArgs(TRUE)
arm <- a[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
say <- function(lab, r) {
  if (inherits(r$fit, "error")) {
    cat(lab, "| ERROR", gsub("\n", " ", conditionMessage(r$fit)), "\n")
    return(invisible(NULL))
  }
  f <- r$fit
  cat(sprintf("%s | ll %.8f conv %d nwarn %d | beta %s | betad %s | theta %s\n",
              lab, as.numeric(logLik(f)), f$opt$convergence, length(r$warn),
              paste(sprintf("%.7g", unlist(f$estimates$beta)), collapse = " "),
              paste(sprintf("%.7g", unlist(f$estimates$betad)), collapse = " "),
              paste(sprintf("%.7g", unlist(f$estimates$theta)), collapse = " ")))
  for (w in r$warn) cat("   W:", substr(gsub("\n", " ", w), 1, 150), "\n")
  # the latent-residual blocks the frame holds, and the mi() rows it keeps
  mp <- f$frame$mi_map
  for (nm in names(mp)) {
    cat("   mi_map", nm, "rows", paste(mp[[nm]]$rows, collapse = ","),
        "unread", paste(mp[[nm]]$unread, collapse = ","), "\n")
  }
}
mk <- function(seed) {
  set.seed(73000 + seed)
  G <- 10L; K <- 6L
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- rnorm(n)
  e <- as.vector(replicate(G, arima.sim(list(ar = 0.6), K, sd = 0.5)))
  eta <- 0.6 + 0.4 * d$x + e
  d$ln <- exp(eta + rnorm(n, 0, 0.3))
  d$ga <- rgamma(n, 5, 5 / exp(eta))
  d$wb <- rweibull(n, 2, exp(eta) / gamma(1.5))
  d$k <- rpois(n, exp(eta))
  d$w2 <- rnorm(n)
  d$xm <- d$x + rnorm(n, 0, 0.2)
  d$zs <- rnorm(n)
  d
}
w0 <- function(d, v) {
  d$w <- ifelse(is.na(d[[v]]), 0, 1)
  d[[v]][is.na(d[[v]])] <- 1
  d
}
fams <- list(ln = lognormal(), ga = Gamma(link = "log"), wb = weibull())
terms <- c(ar = "ar(t, g, cov = TRUE)", ma = "ma(t, g, cov = TRUE)",
           arf = "ar(t, g)")
# 1. unread missing values: mi() against a weight of 0, 3 seeds x 3 families x
#    3 terms
for (seed in 1:3) for (fv in names(fams)) for (tm in names(terms)) {
  d <- mk(seed)
  miss <- sort(sample(nrow(d), 7))
  d[[fv]][miss] <- NA
  f1 <- as.formula(sprintf("%s | mi() ~ x + %s", fv, terms[[tm]]))
  f0 <- as.formula(sprintf("%s | weights(w) ~ x + %s", fv, terms[[tm]]))
  say(sprintf("U %s %s seed %d mi()   ", fv, tm, seed),
      fe_fit(bf(f1), family = fams[[fv]], data = d))
  say(sprintf("U %s %s seed %d weight0", fv, tm, seed),
      fe_fit(bf(f0), family = fams[[fv]], data = w0(d, fv)))
}
d <- mk(9)
miss <- c(3L, 14L, 27L, 33L, 50L)
d$ln[miss] <- NA
# 2. read by another response's mu, and by another response's sigma
say("R mu of w2", fe_fit(bf(ln | mi() ~ x + ar(t, g, cov = TRUE),
                             family = lognormal()) +
                          bf(w2 ~ mi(ln), family = gaussian()) +
                          set_rescor(FALSE), data = d))
say("R sigma of w2", fe_fit(bf(ln | mi() ~ x + ar(t, g, cov = TRUE),
                                family = lognormal()) +
                             bf(w2 ~ 1, sigma ~ mi(ln), family = gaussian()) +
                             set_rescor(FALSE), data = d))
# 3. read by its own sigma
say("R own sigma", fe_fit(bf(ln | mi() ~ x + ar(t, g, cov = TRUE),
                              sigma ~ mi(ln)), family = lognormal(), data = d))
# 4. a reader whose own row is missing where the value is: peeled, so unread
d4 <- d
d4$w2[miss] <- NA
say("P w2 missing on the same rows", fe_fit(
  bf(ln | mi() ~ x + ar(t, g, cov = TRUE), family = lognormal()) +
    bf(w2 | mi() ~ mi(ln), family = gaussian()) + set_rescor(FALSE),
  data = d4))
# 5. mi() on a predictor of the latent-residual response; its own missing
#    values unread. The predictor's imputed values then share the inner
#    problem with the latent residuals
d5 <- d
d5$xm[c(5L, 20L, 44L)] <- NA
say("Q mi(xm) in ln, ln | mi()", fe_fit(
  bf(ln | mi() ~ mi(xm) + ar(t, g, cov = TRUE), family = lognormal()) +
    bf(xm | mi() ~ 1, family = gaussian()) + set_rescor(FALSE), data = d5))
say("Q mi(xm) in ln, no mi() on ln", fe_fit(
  bf(ln ~ mi(xm) + ar(t, g, cov = TRUE), family = lognormal()) +
    bf(xm | mi() ~ 1, family = gaussian()) + set_rescor(FALSE),
  data = transform(d5, ln = exp(0.6 + 0.4 * x + rnorm(60, 0, 0.5)))))
say("Q mi(xm) in sigma of ln", fe_fit(
  bf(ln | mi() ~ x + ar(t, g, cov = TRUE), sigma ~ mi(xm),
     family = lognormal()) +
    bf(xm | mi() ~ 1, family = gaussian()) + set_rescor(FALSE), data = d5))
# 6. poisson: mi() on a count response
dk <- d
dk$k[miss] <- NA
say("K poisson mi()", fe_fit(bf(k | mi() ~ x + ar(t, g, cov = TRUE)),
                             family = poisson(), data = dk))
cat("DONE\n")
