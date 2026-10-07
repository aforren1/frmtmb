# Reviewer of lane optima, claim 2: the probit, cloglog and softit
# log-odds. (a) old against new form inside ONE process (pitfall 21),
# on fits of the existing test fixtures' constructions; (b) value and
# tape derivative at the switch points; (c) RTMB's pnorm(log.p = TRUE)
# derivative in the far tail, against the Mills-ratio series.
#   Rscript dev/optima-rev-links.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
Lnew <- get("frmtmb_links", ns)
Lold <- Lnew
Lold$cloglog$logit_eta <- function(eta) {
  t <- exp(eta)
  log(-expm1(-t)) + t
}
Lold$probit$logit_eta <- function(eta) {
  log(RTMB::pnorm(eta)) - log(RTMB::pnorm(-eta))
}
Lold$softit$logit_eta <- function(eta) log(RTMB::logspace_add(0 * eta, eta))
setL <- function(L) {
  unlockBinding("frmtmb_links", ns)
  assign("frmtmb_links", L, envir = ns)
  lockBinding("frmtmb_links", ns)
}

sim_ord_data <- function(seed = 91, n = 600, K = 4) {
  set.seed(seed)
  x <- rnorm(n)
  eta <- 1.2 * x
  tau <- c(-1, 0.3, 1.4)
  u <- runif(n)
  p <- plogis(outer(rep(1, n), tau) - eta)
  y <- rowSums(u > cbind(p, 1)) + 1L
  data.frame(y = y, x = x)
}
set.seed(17)
db <- data.frame(x = rnorm(250))
db$y <- rbinom(250, 1, pnorm(0.2 + 0.8 * db$x))
dbe <- data.frame(x = rnorm(250))
mu <- pnorm(0.2 + 0.5 * dbe$x)
dbe$y <- rbeta(250, mu * 8, (1 - mu) * 8)
set.seed(4)
dbin <- data.frame(x = rnorm(300), nt = sample(3:12, 300, TRUE))
dbin$y <- rbinom(300, dbin$nt, 1 - exp(-exp(-0.5 + 0.7 * dbin$x)))
set.seed(8)
dh <- sim_ord_data(seed = 8, n = 400)
dh$y[sample(400, 60)] <- 0L
fixtures <- list(
  cum_probit_92 = function() frm(bf(y ~ x) + cumulative(link = "probit"),
                                 data = sim_ord_data(seed = 92)),
  cum_probit = function() frm(bf(y ~ x), family = cumulative("probit"),
                              data = sim_ord_data(n = 300)),
  cum_cloglog = function() frm(bf(y ~ x), family = cumulative("cloglog"),
                               data = sim_ord_data(n = 300)),
  cum_softit = function() frm(bf(y ~ x), family = cumulative("softit"),
                              data = sim_ord_data(n = 300)),
  sratio_probit = function() frm(bf(y ~ x), family = sratio("probit"),
                                 data = sim_ord_data(n = 300)),
  sratio_cloglog = function() frm(bf(y ~ x), family = sratio("cloglog"),
                                  data = sim_ord_data(n = 300)),
  cratio_cloglog = function() frm(bf(y ~ x), family = cratio("cloglog"),
                                  data = sim_ord_data()),
  acat_probit = function() frm(bf(y ~ x), family = acat("probit"),
                               data = sim_ord_data(n = 300)),
  bern_probit_17 = function() frm(bf(y ~ x) + bernoulli("probit"), data = db),
  bern_cloglog = function() frm(bf(y ~ x) + bernoulli("cloglog"), data = db),
  binom_cloglog = function() frm(bf(y | trials(nt) ~ x),
                                 family = binomial("cloglog"), data = dbin),
  beta_probit_17 = function() frm(bf(y ~ x) + Beta("probit"), data = dbe),
  hurdle_cum_probit = function() frm(y ~ x,
                                     family = hurdle_cumulative("probit"),
                                     data = dh),
  cum_probit_cs = function() frm(bf(y ~ cs(x)), family = cumulative("probit"),
                                 data = sim_ord_data(n = 300))
)
ulp <- function(a, b) {
  if (identical(a, b)) return(0)
  abs(a - b) / (.Machine$double.eps * max(abs(a), abs(b)))
}
for (nm in names(fixtures)) {
  setL(Lnew)
  fn <- tryCatch(suppressWarnings(fixtures[[nm]]()), error = function(e) e)
  fn2 <- tryCatch(suppressWarnings(fixtures[[nm]]()), error = function(e) e)
  setL(Lold)
  fo <- tryCatch(suppressWarnings(fixtures[[nm]]()), error = function(e) e)
  fo2 <- tryCatch(suppressWarnings(fixtures[[nm]]()), error = function(e) e)
  setL(Lnew)
  if (!inherits(fn2, "error") && !inherits(fo2, "error")) {
    cat(sprintf("CTL %-18s new twice: logLik %s par %s | old twice: logLik %s par %s\n",
                nm, identical(logLik(fn), logLik(fn2)),
                identical(fn$opt$par, fn2$opt$par),
                identical(logLik(fo), logLik(fo2)),
                identical(fo$opt$par, fo2$opt$par)))
  }
  if (inherits(fn, "error") || inherits(fo, "error")) {
    cat(sprintf("FIX %-18s ERROR new=%s old=%s\n", nm,
                if (inherits(fn, "error")) conditionMessage(fn) else "ok",
                if (inherits(fo, "error")) conditionMessage(fo) else "ok"))
    next
  }
  lln <- as.numeric(logLik(fn))
  llo <- as.numeric(logLik(fo))
  pn <- unname(fn$opt$par)
  po <- unname(fo$opt$par)
  # the two objectives at one point: the old fit's optimum
  vn <- fn$obj$fn(po)
  vo <- fo$obj$fn(po)
  gn <- fn$obj$gr(po)
  go <- fo$obj$gr(po)
  cat(sprintf(paste0("FIX %-18s logLik identical %-5s (diff %.3g) ",
                     "par identical %-5s (max rel %.3g) | at one point: ",
                     "fn ulp %.3g, max|gr diff| %.3g, max|gr| %.3g ",
                     "codes %d/%d evals %s/%s\n"),
              nm, identical(lln, llo), lln - llo, identical(pn, po),
              max(abs(pn - po) / pmax(abs(po), 1e-300)), ulp(vn, vo),
              max(abs(gn - go)), max(abs(go)), fn$opt$convergence,
              fo$opt$convergence, format(fn$opt$evals), format(fo$opt$evals)))
}

# (b) the switch points, on RTMB's tape
cat("\n== switch points (value, tape derivative) ==\n")
sw <- function(lk, x) {
  q <- Lnew[[lk]]$logit_eta
  tp <- RTMB::MakeTape(function(e) q(e), x)
  cbind(x = x, value = tp(x), deriv = diag(tp$jacobian(x)))
}
h <- c(-1e-9, -1e-12, 0, 1e-12, 1e-9)
for (lk in c("cloglog", "softit")) {
  m <- sw(lk, -40 + h)
  cat(lk, "at -40 + h:\n")
  print(format(as.data.frame(m), digits = 17), row.names = FALSE)
  o <- Lold[[lk]]$logit_eta
  tpo <- RTMB::MakeTape(function(e) o(e), -40 + h)
  cat("  old form: values", format(tpo(-40 + h), digits = 17), "\n")
  cat("  old form: derivs", format(diag(tpo$jacobian(-40 + h)),
                                   digits = 17), "\n")
}
for (s in c(1, -1)) {
  m <- sw("probit", s * (1000 + h))
  cat("probit at", s, "* (1000 + h):\n")
  print(format(as.data.frame(m), digits = 17), row.names = FALSE)
}
# the exact value and derivative of the probit's log-odds for |x| >= 30
mills <- function(x) {
  a <- abs(x)
  # log Phi(-a) and phi(a) / Phi(-a) from the asymptotic series
  s <- 1 - 1 / a^2 + 3 / a^4 - 15 / a^6 + 105 / a^8
  lphi_neg <- -a^2 / 2 - log(a) - 0.5 * log(2 * pi) + log(s)
  dm <- a / s
  list(lpn = lphi_neg, d = dm)
}
cat("\nprobit derivative jump at |eta| = 1000: left (exact)",
    format(mills(1000)$d, digits = 17), "; right (continuation) 1000",
    "; relative", format((mills(1000)$d - 1000) / 1000, digits = 3), "\n")

# (c) RTMB's pnorm(log.p = TRUE) derivative, far tail
cat("\n== RTMB pnorm(x, log.p = TRUE): value and derivative at x < 0 ==\n")
xs <- -c(30, 100, 1e3, 1e4, 3e4, 1e5, 1e6, 3e6, 9.4e6, 1e7, 1e8, 1e9,
         4.46e9, 1e10, 1e100, 1e154, 1.4e154)
tp <- RTMB::MakeTape(function(x) RTMB::pnorm(x, log.p = TRUE), xs)
v <- tp(xs)
d <- diag(tp$jacobian(xs))
for (i in seq_along(xs)) {
  m <- mills(xs[i])
  cat(sprintf("x %10.4g  value rel err %10.3g  deriv rel err %10.3g\n",
              xs[i], (v[i] - m$lpn) / abs(m$lpn), (d[i] - m$d) / m$d))
}
# a dense scan near 9.4e6 and 4.46e9
for (c0 in c(9.4e6, 4.46e9)) {
  x <- -seq(c0 * 0.999, c0 * 1.001, length.out = 2001)
  tp <- RTMB::MakeTape(function(x) RTMB::pnorm(x, log.p = TRUE), x)
  d <- diag(tp$jacobian(x))
  ex <- vapply(x, function(z) mills(z)$d, 0)
  re <- abs(d - ex) / ex
  cat(sprintf("scan near %.3g: non-finite %d of 2001; max rel err %.3g; median %.3g\n",
              c0, sum(!is.finite(d)), max(re[is.finite(re)]),
              stats::median(re[is.finite(re)])))
}
# the same scans through the lane's probit log-odds (clamped at 1000)
q <- Lnew$probit$logit_eta
x <- -c(seq(1001, 1e4, length.out = 500), seq(9.39e6, 9.41e6, length.out = 500),
        seq(4.455e9, 4.465e9, length.out = 1001))
tp <- RTMB::MakeTape(function(e) q(e), x)
d <- diag(tp$jacobian(x))
cat("lane probit log-odds far tail: non-finite derivative", sum(!is.finite(d)),
    "of", length(x), "; non-finite value", sum(!is.finite(tp(x))), "\n")
