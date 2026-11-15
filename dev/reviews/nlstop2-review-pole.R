# Review of lane nlstop2, attack 1: the pole escape's scope, triggers,
# and bit-identity of fits clear of poles; attack 2: the limit check on
# the explin ridge at either ending.
# REV_ARM=before|after Rscript dev/reviews/nlstop2-review-pole.R
source("dev/reviews/nlstop2-review-common.R")
mm_data <- function(seed, family) {
  set.seed(seed * 1000 + 1)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                  xq = round(runif(n, 0, 5), 3),
                  xb = round(runif(n, 0, 3000), 1),
                  g = factor(sample(1:15, n, TRUE)))
  d$ug <- rnorm(15, 0, 0.4)[d$g]
  lat <- 3 * d$xp / (1 + d$xp)
  d$nt <- 10L
  d$y <- switch(family,
                bernoulli = rbinom(n, 1, plogis(lat - 2)),
                binomial = rbinom(n, 10, plogis(lat - 2)),
                gaussian = lat + rnorm(n, 0, 0.5))
  d
}
cnt <- new.env()
# count objective evaluations of every fit by tracing optimize_obj? Use
# the time and the message classes instead; cost is the lane's claim
row <- function(tag, expr) {
  t0 <- proc.time()[[3]]
  r <- tt(expr)
  el <- proc.time()[[3]] - t0
  if (inherits(r$v, "tterr")) {
    cat("ROW", tag, "ERROR", substr(r$v$msg, 1, 200), "\n")
    return(invisible(NULL))
  }
  f <- r$v
  b <- f$estimates$beta
  cat("ROW", tag, "ll", format(llv(f), digits = 17), "code",
      f$opt$convergence, "beta", paste(format(b, digits = 10),
                                        collapse = " "),
      "warn", if (length(r$w)) paste(unique(r$w), collapse = "+") else "-",
      "msg", if (length(r$m)) paste(unique(r$m), collapse = "+") else "-",
      "t", round(el, 2), "\n")
  se <- tryCatch(sqrt(diag(suppressWarnings(vcov(f)))), error = function(e) NA)
  cat("   se", paste(format(se, digits = 10), collapse = " "), "\n")
  invisible(f)
}

# 1. the lane's own rows, and clear rows for bit identity
for (s in c(11, 44)) {
  d <- mm_data(s, "binomial")
  row(paste("mm binomial", s), frm(bf(y | trials(nt) ~ Vm * xp / (K + xp),
                                      Vm + K ~ 1, nl = TRUE),
                                   family = binomial(), data = d))
}
for (s in c(11, 22, 55)) {
  d <- mm_data(s, "gaussian")
  row(paste("mm gaussian", s), frm(bf(y ~ Vm * xp / (K + xp), Vm + K ~ 1,
                                      nl = TRUE), data = d))
}
# 2. K with a group-level term
for (s in c(11, 44)) {
  d <- mm_data(s, "binomial")
  row(paste("mm binomial K~(1|g)", s),
      frm(bf(y | trials(nt) ~ Vm * xp / (K + xp), Vm ~ 1, K ~ 1 + (1 | g),
             nl = TRUE), family = binomial(), data = d))
}
# 3. coefficient enters the denominator nonlinearly: K^2 + xp, exp(K) + xp
for (s in c(11, 44)) {
  d <- mm_data(s, "binomial")
  row(paste("mm binomial K^2", s),
      frm(bf(y | trials(nt) ~ Vm * xp / (K^2 + xp), Vm + K ~ 1, nl = TRUE),
          family = binomial(), data = d))
}
# 4. two coefficients in one denominator: a / (b + c * xp)
for (s in c(11, 44)) {
  d <- mm_data(s, "binomial")
  row(paste("two coef a/(b+c*xp)", s),
      frm(bf(y | trials(nt) ~ a * xp / (b + c * xp), a + b + c ~ 1,
             nl = TRUE), family = binomial(), data = d))
}
# 5. pole inside exp()
for (s in c(11, 44)) {
  d <- mm_data(s, "binomial")
  row(paste("exp pole", s),
      frm(bf(y | trials(nt) ~ Vm * exp(-1 / (K + xp)), Vm + K ~ 1,
             nl = TRUE), family = binomial(), data = d))
}
# 6. negative literal power
for (s in c(11, 44)) {
  d <- mm_data(s, "binomial")
  row(paste("neg power", s),
      frm(bf(y | trials(nt) ~ Vm * xp * (K + xp)^-1, Vm + K ~ 1,
             nl = TRUE), family = binomial(), data = d))
}
# 7. a/(1+b*x) with b near 1/max|x|: walls close to the fit, all
# denominators positive, true model
for (s in 1:5) {
  set.seed(s)
  d <- data.frame(x = rnorm(300))
  bt <- 0.95 / max(abs(d$x))
  d$y <- 2 / (1 + bt * d$x) + rnorm(300, 0, 0.3)
  row(paste("a/(1+bx) near wall", s, "btrue", signif(bt, 3)),
      frm(bf(y ~ a / (1 + b * x), a + b ~ 1, nl = TRUE), data = d))
}
# 8. true pole inside the data range: y = 1/(x + b), b = -2, x in [0, 4]
# away from 2; the maximum is between walls (mixed signs)
for (s in 1:3) {
  set.seed(s)
  x <- runif(400, 0, 4)
  x <- x[abs(x - 2) > 0.3]
  d <- data.frame(x = x)
  d$y <- 1 / (d$x - 2) + rnorm(nrow(d), 0, 0.2)
  row(paste("true pole inside", s),
      frm(bf(y ~ a / (x + b), a + b ~ 1, nl = TRUE), data = d,
          start = list(beta = c(1, -1.9))))
  row(paste("true pole inside default", s),
      frm(bf(y ~ a / (x + b), a + b ~ 1, nl = TRUE), data = d))
}
# 9. true denominator never changes sign but sits near a wall
for (s in 1:3) {
  set.seed(s)
  d <- data.frame(x = runif(200, 0.1, 5))
  d$y <- 0.5 / (d$x - 0.095) + rnorm(200, 0, 0.3)
  row(paste("near wall true", s),
      frm(bf(y ~ a / (x + b), a + b ~ 1, nl = TRUE), data = d))
}
# 10. two nlpars bodies (dpar sigma nonlinear too)
for (s in c(11, 44)) {
  d <- mm_data(s, "gaussian")
  row(paste("nl mu+sigma", s),
      frm(bf(y ~ Vm * xp / (K + xp), Vm + K ~ 1,
             nlf(sigma ~ s0 + s1 / (s2 + xp)), s0 + s1 + s2 ~ 1, nl = TRUE),
          data = d))
}
# 11. restarts = 0, custom optimizer, given start, sign probe off
d <- mm_data(11, "binomial")
fm <- bf(y | trials(nt) ~ Vm * xp / (K + xp), Vm + K ~ 1, nl = TRUE)
row("mm11 restarts=0", frm(fm, family = binomial(), data = d,
                           control = frmtmb_control(restarts = 0)))
row("mm11 optim", frm(fm, family = binomial(), data = d,
                      control = frmtmb_control(optimizer = "optim")))
options(frmtmb.nl_sign_probe = FALSE)
row("mm11 probe off", frm(fm, family = binomial(), data = d))
options(frmtmb.nl_sign_probe = TRUE)
row("mm11 start K -2.5", frm(fm, family = binomial(), data = d,
                             start = list(beta = c(-0.1, -2.5))))
# 12. explin ridge: the limit check at either ending, over ulp draws
set.seed(2)
d <- data.frame(x = runif(1000, 0, 3))
d$y <- 2 * d$x + rnorm(1000, 0, 0.5)
for (k in 0:11) {
  dk <- d
  dk$y <- dk$y * (1 + 2e-15 * k)
  r <- tt(frm(bf(y ~ a * (exp(b * x) - 1), a + b ~ 1, nl = TRUE), data = dk,
              start = list(beta = c(2, 0.5))))
  f <- r$v
  msg <- if (!inherits(f, "tterr")) tryCatch(ns$nl_limit_message(f),
                                             error = function(e) "ERR")
  cat("ROW explin k", k, "code", f$opt$convergence, "ll",
      format(llv(f), digits = 12), "warn",
      paste(unique(r$w), collapse = "+"), "direct_check",
      if (is.null(msg)) "NULL" else substr(msg, 1, 60), "\n")
}
cat("DONE\n")
