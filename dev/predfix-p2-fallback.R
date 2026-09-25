# Punch round 2, requirement (c): on the calibration sets, does any fit
# that engaged before now fall back and lose its improvement? Refits the
# DEFAULT on the engaged cells of the 945-fit scale scan and of the
# random-slope calibration (parts 1 and 2), and compares with the
# FALSE and TRUE values those runs saved.
#   PREDFIX_ARM=lane Rscript dev/predfix-p2-fallback.R > dev/predfix-log/p2-fallback.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
dl <- "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/"
q <- function(expr) suppressWarnings(suppressMessages(expr))
fitd <- function(...) {
  f <- try(q(frm(...)), silent = TRUE)
  if (inherits(f, "try-error")) list(ll = NA, engaged = NA) else
    list(ll = as.numeric(logLik(f)), engaged = !is.null(f$par_units))
}

## the scale scan: scales below 1e-3 engage
scan <- readRDS(paste0(dl, "scalescan.rds"))
mk <- function(seed, n = 250) {
  set.seed(seed)
  x <- rnorm(n)
  z <- runif(n)
  data.frame(
    x = x, z = z, f = factor(sample(letters[1:3], n, TRUE)),
    pois = rpois(n, exp(0.5 + 0.4 * x)),
    bin = rbinom(n, 1, plogis(0.3 + 0.8 * x)),
    gau = 1 + 0.5 * x + rnorm(n),
    gam = rgamma(n, 2, 2 / exp(0.2 + 0.3 * x)),
    nb = rnbinom(n, mu = exp(0.5 + 0.4 * x), size = 3))
}
famof <- list(pois = poisson(), bin = bernoulli(), gau = gaussian(),
              gam = Gamma(link = "log"), nb = negbinomial())
eng <- scan[scan$scale < 1e-3, ]
eng$new_default <- NA_real_
eng$engaged <- NA
for (seed in unique(eng$seed)) {
  d <- mk(seed)
  for (i in which(eng$seed == seed)) {
    d$xs <- d$x * eng$scale[i]
    r <- fitd(bf(stats::as.formula(paste(eng$y[i], "~", eng$rhs[i]))),
              family = famof[[eng$y[i]]], data = d)
    eng$new_default[i] <- r$ll
    eng$engaged[i] <- r$engaged
  }
}
eng$lost <- eng$auto - eng$new_default
cat("scale scan, engaged cells:", nrow(eng), "; default engaged:",
    sum(eng$engaged, na.rm = TRUE), "; fell back:",
    sum(!eng$engaged, na.rm = TRUE), "; errored:", sum(is.na(eng$engaged)),
    "\n  below autoscale = TRUE by more than 1e-6:",
    sum(eng$lost > 1e-6, na.rm = TRUE), "; largest loss",
    signif(max(eng$lost, na.rm = TRUE), 4), "\n")

## the random-slope calibration, parts 1 and 2 (spreads below 0.05)
mk2 <- function(seed, tau, fam) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x <- rnorm(300)
  eta <- 0.5 + 0.4 * x + rnorm(20, 0, 0.6)[g] + rnorm(20, 0, tau)[g] * x
  y <- if (fam == "gau") eta + rnorm(300) else rpois(300, exp(eta))
  data.frame(g = g, x0 = x, y = y)
}
sl <- rbind(readRDS(paste0(dl, "p1-slopecal.rds")),
            readRDS(paste0(dl, "p1-slopecal2.rds")))
sl <- sl[sl$s < 0.05, ]
sl$short_default <- NA_real_
sl$engaged <- NA
fo <- bf(y ~ x + (1 + x | g))
for (k in seq_len(nrow(sl))) {
  d <- mk2(sl$seed[k], sl$tau[k], sl$fam[k])
  d$x <- d$x0 * sl$s[k]
  r <- fitd(fo, family = if (sl$fam[k] == "gau") gaussian() else poisson(),
            data = d)
  # the reference is the saved one: ref = ll_FALSE + short_false, and
  # the saved short_true is TRUE's shortfall against it
  ll_false <- NA
  sl$engaged[k] <- r$engaged
  sl$new_ll[k] <- r$ll
}
# recover ref from a fresh FALSE fit is costly; the saved shortfalls are
# against ref, so compare the default to TRUE through them
sl$ref_minus_default <- NA_real_
for (k in seq_len(nrow(sl))) {
  d <- mk2(sl$seed[k], sl$tau[k], sl$fam[k])
  d$x <- d$x0 * sl$s[k]
  lf <- fitd(fo, family = if (sl$fam[k] == "gau") gaussian() else poisson(),
             data = d, control = frmtmb_control(autoscale = FALSE))$ll
  ref <- lf + sl$short_false[k]
  sl$ref_minus_default[k] <- ref - sl$new_ll[k]
}
sl$lost <- sl$ref_minus_default - sl$short_true
cat("slope calibration, spread below 0.05:", nrow(sl), "fits; default engaged:",
    sum(sl$engaged, na.rm = TRUE), "; fell back:", sum(!sl$engaged, na.rm = TRUE),
    "; errored:", sum(is.na(sl$engaged)),
    "\n  below autoscale = TRUE by more than 1e-6:",
    sum(sl$lost > 1e-6, na.rm = TRUE), "; largest loss",
    signif(max(sl$lost, na.rm = TRUE), 4),
    "\n  short of the reference by more than 1e-3:",
    sum(sl$ref_minus_default > 1e-3, na.rm = TRUE), "\n")
if (any(sl$lost > 1e-6, na.rm = TRUE)) print(sl[which(sl$lost > 1e-6), ])
saveRDS(list(scan = eng, slope = sl), paste0(dl, "p2-fallback.rds"))
