# Area 6: the 20 Stan functions on the tape: value and derivative at extremes,
# against an exact mpfr-free reference (log-space formulas evaluated by hand) and
# against Stan where it can be run. Seed 8601.
source("dev/review-nlfix-common.R")
ns <- asNamespace("frmtmb")
F <- ns$nl_stan_funs
grad_tape <- function(f, x) {
  # value and gradient of f on the tape, by RTMB's own derivative
  tp <- RTMB::MakeTape(function(p) {
    do.call(f, as.list(p))
  }, x)
  list(v = as.numeric(tp(x)), g = as.numeric(tp$jacobian(x)))
}
fd <- function(f, x, i) {
  h <- 1e-6 * max(1, abs(x[i]))
  xp <- x; xm <- x
  xp[i] <- x[i] + h; xm[i] <- x[i] - h
  (do.call(f, as.list(xp)) - do.call(f, as.list(xm))) / (2 * h)
}
ref <- list(
  logit = list(v = function(x) qlogis(x), g = function(x) 1 / (x * (1 - x))),
  inv_logit = list(v = function(x) plogis(x), g = function(x) dlogis(x)),
  inv_cloglog = list(v = function(x) -expm1(-exp(x)), g = function(x) exp(x - exp(x))),
  Phi = list(v = function(x) pnorm(x), g = function(x) dnorm(x)),
  square = list(v = function(x) x^2, g = function(x) 2 * x),
  cbrt = list(v = function(x) sign(x) * abs(x)^(1/3), g = function(x) 1 / (3 * abs(x)^(2/3))),
  exp2 = list(v = function(x) 2^x, g = function(x) log(2) * 2^x),
  inv = list(v = function(x) 1 / x, g = function(x) -1 / x^2),
  inv_sqrt = list(v = function(x) 1 / sqrt(x), g = function(x) -0.5 * x^-1.5),
  inv_square = list(v = function(x) 1 / x^2, g = function(x) -2 / x^3),
  log1m = list(v = function(x) log1p(-x), g = function(x) -1 / (1 - x)),
  step = list(v = function(x) as.numeric(x >= 0), g = function(x) 0),
  log1p_exp = list(v = function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x))), g = function(x) plogis(x)),
  log1m_exp = list(v = function(x) ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x))), g = function(x) -1 / expm1(-x)),
  log_inv_logit = list(v = function(x) plogis(x, log.p = TRUE), g = function(x) plogis(-x)),
  log1m_inv_logit = list(v = function(x) plogis(-x, log.p = TRUE), g = function(x) -plogis(x))
)
grids <- list(
  logit = c(1e-300, 1e-20, 1e-8, 0.5, 1 - 1e-8, 1 - 1e-15),
  inv_logit = c(-800, -745, -710, -40, -1, 0, 1, 40, 710, 745, 800, 1e4),
  inv_cloglog = c(-800, -40, -1, 0, 1, 5, 6, 7, 10, 500, 710, 800),
  Phi = c(-40, -38, -37, -9, 0, 9, 37, 38, 40),
  square = c(-1e160, 1e-170, 0, 1e160),
  cbrt = c(-8, -1e-300, 0, 1e-300, 8, 1e300),
  exp2 = c(-1100, -1, 0, 1, 1023),
  inv = c(-1e-300, 1e-300, 1e300, 0),
  inv_sqrt = c(1e-300, 1, 1e300, 0),
  inv_square = c(1e-160, 1, 1e160, 0),
  log1m = c(-1e300, -1, 1e-20, 0.5, 1 - 1e-15),
  step = c(-1e300, -1e-300, 0, 1e-300, 1e300),
  log1p_exp = c(-800, -745, -40, -1, 0, 1, 40, 710, 800, 1e4),
  log1m_exp = c(-800, -40, -1, -1e-5, -1e-20, -1e-300),
  log_inv_logit = c(-1e4, -800, -40, 0, 40, 800, 1e4),
  log1m_inv_logit = c(-1e4, -800, -40, 0, 40, 800, 1e4)
)
cat(sprintf("%-16s %10s %14s %14s %14s %14s %s\n", "fun", "x", "dbl", "tape", "ref value", "tape grad", "ref grad"))
bad <- 0L
for (nm in names(grids)) {
  for (x in grids[[nm]]) {
    dbl <- tryCatch(suppressWarnings(as.numeric(F[[nm]](x))), error = function(e) NA_real_)
    tp <- tryCatch(suppressWarnings(grad_tape(F[[nm]], x)), error = function(e) list(v = NA, g = NA))
    rv <- suppressWarnings(ref[[nm]]$v(x))
    rg <- suppressWarnings(ref[[nm]]$g(x))
    flag <- character(0)
    rel <- function(a, b) if (!is.finite(a) || !is.finite(b)) { if (identical(a, b) || (is.nan(a) && is.nan(b))) 0 else Inf } else abs(a - b) / max(1e-300, abs(b))
    if (rel(dbl, rv) > 1e-10) flag <- c(flag, "VAL-dbl")
    if (rel(tp$v, rv) > 1e-10) flag <- c(flag, "VAL-tape")
    if (rel(tp$g, rg) > 1e-8) flag <- c(flag, "GRAD")
    if (length(flag)) bad <- bad + 1L
    cat(sprintf("%-16s %10.4g %14.7g %14.7g %14.7g %14.7g %14.7g %s\n", nm, x, dbl, tp$v, rv, tp$g, rg, paste(flag, collapse = ",")))
  }
}
cat("two-argument functions\n")
two <- list(
  log_sum_exp = list(g = rbind(c(800, 800), c(-800, -800), c(800, -800), c(-800, 800), c(0, 0), c(1e4, 1e4), c(710, 709)),
    v = function(x, y) pmax(x, y) + log1p(exp(pmin(x, y) - pmax(x, y))),
    gx = function(x, y) 1 / (1 + exp(y - x)), gy = function(x, y) 1 / (1 + exp(x - y))),
  log_diff_exp = list(g = rbind(c(800, 799), c(-800, -801), c(1e-3, 0), c(1e-12, 0), c(0, -800), c(800, -800), c(1, 1)),
    v = function(x, y) x + ifelse(y - x > -log(2), log(-expm1(y - x)), log1p(-exp(y - x))),
    gx = function(x, y) -1 / expm1(y - x), gy = function(x, y) 1 / expm1(y - x)),
  hypot = list(g = rbind(c(3, 4), c(1e200, 1e200), c(1e-200, 1e-200), c(0, 0), c(0, 3)),
    v = function(x, y) { m <- pmax(abs(x), abs(y)); m * sqrt(1 + (pmin(abs(x), abs(y)) / m)^2) },
    gx = function(x, y) { m <- pmax(abs(x), abs(y)); (x / m) / sqrt(1 + (pmin(abs(x), abs(y)) / m)^2) },
    gy = function(x, y) { m <- pmax(abs(x), abs(y)); (y / m) / sqrt(1 + (pmin(abs(x), abs(y)) / m)^2) }),
  pow = list(g = rbind(c(2, 3), c(0, 2), c(0, 0), c(0, 0.5), c(4, -0.5), c(-2, 2)),
    v = function(x, y) x^y, gx = function(x, y) y * x^(y - 1), gy = function(x, y) if (x == 0) 0 else x^y * log(x))
)
for (nm in names(two)) {
  G <- two[[nm]]$g
  for (i in seq_len(nrow(G))) {
    x <- G[i, ]
    dbl <- tryCatch(suppressWarnings(as.numeric(F[[nm]](x[1], x[2]))), error = function(e) NA_real_)
    tp <- tryCatch(suppressWarnings(grad_tape(F[[nm]], x)), error = function(e) list(v = NA, g = c(NA, NA)))
    cat(sprintf("%-14s (%8.4g,%8.4g) dbl %12.7g tape %12.7g ref %12.7g | grad tape %12.7g %12.7g ref %12.7g %12.7g\n", nm, x[1], x[2], dbl, tp$v,
                two[[nm]]$v(x[1], x[2]), tp$g[1], tp$g[2], two[[nm]]$gx(x[1], x[2]), two[[nm]]$gy(x[1], x[2])))
  }
}
cat("\nflagged rows (one-arg):", bad, "\n")
cat("\n--- a fit that reaches the extreme: saturating bodies\n")
set.seed(8601)
n <- 300
d <- data.frame(x = rnorm(n), z = runif(n, 1, 5))
d$yb <- rbinom(n, 1, plogis(40 * d$x))
d$yc <- rbinom(n, 1, 1 - exp(-exp(3 * d$x)))
d$yg <- plogis(30 * d$x) + rnorm(n, 0, 0.2)
d$yl <- rnorm(n, log1p(exp(3 * d$x)), 0.5)
fit1 <- function(lab, f, fam, d, ...) {
  r <- tt(frm(f, family = fam, data = d, ...))
  if (iserr(r)) cat(sprintf("%-28s ERROR %s\n", lab, substr(r$msg, 1, 200)))
  else cat(sprintf("%-28s ll=%.5f conv=%s beta=%s WARN[%s]\n", lab, as.numeric(logLik(r)), r$opt$convergence,
                   paste(format(r$estimates$beta, digits = 5), collapse = ","), if (length(attr(r, "w"))) substr(attr(r, "w")[1], 1, 140) else ""))
  invisible(r)
}
fit1("inv_logit separation", bf(yb ~ inv_logit(a * x), a ~ 1, nl = TRUE), bernoulli("identity"), d, start = list(beta = c("a_(Intercept)" = 1)))
fit1("plogis linear, same", bf(yb ~ a * x, a ~ 1, nl = TRUE), bernoulli(), d)
fit1("inv_cloglog body", bf(yc ~ inv_cloglog(a + b * x), a + b ~ 1, nl = TRUE), bernoulli("identity"), d, start = list(beta = c("a_(Intercept)" = 0.1, "b_(Intercept)" = 1)))
fit1("cloglog linear", yc ~ x, bernoulli("cloglog"), d)
fit1("gauss inv_logit steep", bf(yg ~ inv_logit(a * x), a ~ 1, nl = TRUE), gaussian(), d, start = list(beta = c("a_(Intercept)" = 1)))
fit1("log1p_exp", bf(yl ~ log1p_exp(a * x), a ~ 1, nl = TRUE), gaussian(), d, start = list(beta = c("a_(Intercept)" = 1)))
fit1("log1p_exp start 300", bf(yl ~ log1p_exp(a * x), a ~ 1, nl = TRUE), gaussian(), d, start = list(beta = c("a_(Intercept)" = 300)))
fit1("inv_logit start 300", bf(yg ~ inv_logit(a * x), a ~ 1, nl = TRUE), gaussian(), d, start = list(beta = c("a_(Intercept)" = 300)))
fit1("inv_cloglog start 300", bf(yg ~ inv_cloglog(a * x), a ~ 1, nl = TRUE), gaussian(), d, start = list(beta = c("a_(Intercept)" = 300)))
fit1("Phi steep start 300", bf(yg ~ Phi(a * x), a ~ 1, nl = TRUE), gaussian(), d, start = list(beta = c("a_(Intercept)" = 300)))
fit1("cbrt default start", bf(yl ~ cbrt(a * x), a ~ 1, nl = TRUE), gaussian(), d)
fit1("cbrt start 1", bf(yl ~ cbrt(a * x), a ~ 1, nl = TRUE), gaussian(), d, start = list(beta = c("a_(Intercept)" = 1)))
fit1("step, change point", bf(yg ~ a * step(x - c), a + c ~ 1, nl = TRUE), gaussian(), d, start = list(beta = c("a_(Intercept)" = 1, "c_(Intercept)" = 0.1)))
cat("\n--- x^b with a zero row in the data: gaussian, before and after\n")
d$xz <- d$z; d$xz[1:5] <- 0
d$yp <- 2 * d$xz^0.5 + rnorm(n, 0, 0.3)
fit1("a * xz^b", bf(yp ~ a * xz^b, a + b ~ 1, nl = TRUE), gaussian(), d, start = list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 1)))
fit1("a * pow(xz, b)", bf(yp ~ a * pow(xz, b), a + b ~ 1, nl = TRUE), gaussian(), d, start = list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 1)))
keep <- d$xz > 0
fit1("a * z^b (no zero rows)", bf(yp ~ a * xz^b, a + b ~ 1, nl = TRUE), gaussian(), d[keep, ], start = list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 1)))
cat("\n--- names outside the 20\n")
for (fn in c("cloglog", "logm1", "expp1", "logit_scaled", "inv_logit_scaled", "multiply_log", "log_mean_exp", "log_expm1", "scale_unit", "fabs", "inv_odds", "softit", "inv_softit", "log_inv_softit", "incgamma", "softmax", "log_softmax", "int_step", "lbeta", "Phi_approx", "fmin", "fmax", "lchoose", "lgamma", "digamma", "log_mix", "binomial_logit_lpmf", "abs", "sign", "pnorm", "expm1", "log1p", "atan2", "sinh")) {
  body <- if (fn %in% c("multiply_log", "lbeta", "fmin", "fmax", "lchoose", "log_mix", "atan2")) sprintf("yg ~ %s(a, x + 3)", fn) else sprintf("yg ~ %s(a + x + 3)", fn)
  r <- tt(frm(bf(as.formula(body), a ~ 1, nl = TRUE), gaussian(), data = d, start = list(beta = c("a_(Intercept)" = 0.2))))
  cat(sprintf("%-18s %s\n", fn, if (iserr(r)) paste("ERR", substr(r$msg, 1, 150)) else paste0("fits ll=", format(as.numeric(logLik(r)), digits = 7))))
  sc <- tryCatch(brms::stancode(brms::bf(as.formula(body), a ~ 1, nl = TRUE), data = d), error = function(e) NULL)
}
cat("\n--- shadow rule\n")
d$inv_logit <- d$x
d$step <- d$x
r <- tt(frm(bf(yg ~ a + inv_logit, a ~ 1, nl = TRUE), gaussian(), data = d)); show("column inv_logit read as data", if (iserr(r)) r else r$estimates$beta)
r <- tt(frm(bf(yg ~ inv_logit(a * x), a ~ 1, nl = TRUE), gaussian(), data = d, start = list(beta = c("a_(Intercept)" = 1)))); show("column and call", if (iserr(r)) r else r$estimates$beta)
r <- tt(frm(bf(yg ~ inv_logit(a * inv_logit), a ~ 1, nl = TRUE), gaussian(), data = d, start = list(beta = c("a_(Intercept)" = 1)))); show("column AND a call to the same name in one body", if (iserr(r)) r else r$estimates$beta)
d$inv_logit <- NULL; d$step <- NULL
f1 <- local({ inv_logit <- function(x) x * 0 + 0.5 + 0.1 * x; bf(yg ~ inv_logit(a * x), a ~ 1, nl = TRUE) })
r <- tt(frm(f1, gaussian(), data = d, start = list(beta = c("a_(Intercept)" = 1)))); show("user function in formula env", if (iserr(r)) r else r$estimates$beta)
attach(list(logit = function(x) 2 * x), name = "rvattached", warn.conflicts = FALSE)
d$p <- plogis(d$x)
d$yq <- 0.5 * qlogis(d$p) + rnorm(n, 0, 0.2)
r <- tt(frm(bf(yq ~ a * logit(p), a ~ 1, nl = TRUE), gaussian(), data = d)); show("logit from an attached package (2 x)", if (iserr(r)) r else r$estimates$beta)
detach("rvattached")
r <- tt(frm(bf(yq ~ a * logit(p), a ~ 1, nl = TRUE), gaussian(), data = d)); show("logit, Stan's", if (iserr(r)) r else r$estimates$beta)
cat("\n--- names do not leak\n")
r <- tt(frm(yg ~ inv_logit(x), data = d)); show("linear formula", if (iserr(r)) r else 1)
r <- tt(frm(bf(yg ~ a * x, a ~ 1 + inv_logit(z), nl = TRUE), data = d)); show("nlpar formula (linear)", if (iserr(r)) r else 1)
r <- tt(frm(bf(yg ~ inv_logit(x), sigma ~ inv_logit(z)), data = d)); show("dpar formula", if (iserr(r)) r else 1)
print(exists("inv_logit", envir = globalenv()))
print("inv_logit" %in% getNamespaceExports("frmtmb"))
fb <- tt(brms::stancode(brms::bf(yg ~ inv_logit(x)), data = d)); show("brms linear formula", if (iserr(fb)) fb else 1)
