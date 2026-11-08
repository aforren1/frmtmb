# Review of lane nlstart: nl_par_units(), against the linear spelling.
# Usage: REV_ARM=after|before Rscript dev/review-nlstart-units.R
source("dev/review-nlstart-common.R")
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
show <- function(label, fb, fl, cb, cl) {
  eb <- unname(fb$estimates$beta[cb]); el <- unname(fl$estimates$beta[cl])
  vb <- sqrt(diag(vcov(fb))); vl <- sqrt(diag(vcov(fl)))
  cat("  vcov names body:", names(vb), "; lin:", names(vl), "\n")
  sb <- unname(vb[seq_along(cb)]); sl <- unname(vl[seq_along(cl)])
  cat(sprintf("%s: logLik body %.8f lin %.8f | est rel %.2e | se rel %.2e\n",
              label, llv(fb), llv(fl), rel(eb, el), rel(sb, sl)))
  wb <- tryCatch(confint(fb), error = function(e) conditionMessage(e))
  wl <- tryCatch(confint(fl), error = function(e) conditionMessage(e))
  if (is.matrix(wb) && is.matrix(wl)) {
    cat(sprintf("  Wald confint rel %.2e\n",
                rel(unname(wb[cb, , drop = FALSE]), unname(wl[cl, , drop = FALSE]))))
  } else cat("  Wald confint:", substr(paste(wb), 1, 120), "|",
             substr(paste(wl), 1, 120), "\n")
  pb <- tryCatch(confint(fb, parm = cb, method = "profile"),
                 error = function(e) conditionMessage(e))
  pl <- tryCatch(confint(fl, parm = cl, method = "profile"),
                 error = function(e) conditionMessage(e))
  if (is.matrix(pb) && is.matrix(pl)) {
    cat(sprintf("  profile confint rel %.2e  body %s  lin %s\n",
                rel(unname(pb[, 1:2]), unname(pl[, 1:2])),
                paste(format(pb[, 1:2], digits = 6), collapse = " "),
                paste(format(pl[, 1:2], digits = 6), collapse = " ")))
  } else cat("  profile confint:", substr(paste(pb), 1, 150), "|",
             substr(paste(pl), 1, 150), "\n")
  cat(sprintf("  fitted rel %.2e  frm_linpred rel %.2e\n",
              rel(fitted(fb)[, 1], fitted(fl)[, 1]),
              rel(as.numeric(frm_linpred(fb)), as.numeric(frm_linpred(fl)))))
  sm <- tryCatch(summary(fb)$coefficients, error = function(e) NULL)
  if (!is.null(sm)) {
    cat("  summary body rows:\n")
    print(sm)
  }
  pt <- tryCatch(par_template(fb), error = function(e) conditionMessage(e))
  cat("  par_template(fit) beta:", paste(format(unlist(pt$beta), digits = 6),
                                         collapse = " "), "\n")
  cat("  par_units:", paste(format(fb$par_units %||% 1, digits = 4),
                            collapse = " "), "\n")
}

cat("=== tiny x (seed 6201)\n")
set.seed(6201)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- 0.5 * d$x + rnorm(n)
d$xt <- d$x * 1e-20
fb <- tt(frm(bf(y ~ a * xt, a ~ 1, nl = TRUE), data = d))
fl <- tt(frm(y ~ 0 + xt, data = d))
show("a * xt vs y ~ 0 + xt", fb, fl, "a_(Intercept)", "xt")
cat("  warnings body:", attr(fb, "wc"), "\n")

cat("=== prior normal(0, 1e-3) on a, tiny x\n")
pb <- tt(frm(bf(y ~ a * xt, a ~ 1, nl = TRUE), data = d,
             prior = prior(normal(0, 1e-3), nlpar = a)))
pl <- tt(frm(y ~ 0 + xt, data = d, prior = prior(normal(0, 1e-3), class = b)))
cat(sprintf("body a %.6g  lin b %.6g  logLik body %.8f lin %.8f\n",
            pb$estimates$beta[[1]], pl$estimates$beta[[1]], llv(pb), llv(pl)))
# the exact MAP of a given the fit's sigma: (x'y / s^2) / (x'x / s^2 + 1/tau^2)
s <- exp(pb$estimates$betad[[1]])
amap <- (sum(d$xt * d$y) / s^2) / (sum(d$xt^2) / s^2 + 1 / 1e-3^2)
cat(sprintf("exact MAP of a at the fit's sigma: %.6g\n", amap))
cat(sprintf("sigma body %.10f lin %.10f; conv body %s lin %s; iter body %s lin %s\n",
            pb$estimates$betad[[1]], pl$estimates$betad[[1]],
            pb$opt$convergence, pl$opt$convergence, pb$opt$iterations,
            pl$opt$iterations))
cat("  body units:", format(pb$par_units %||% 1, digits = 4),
    " grad at opt:", format(pb$obj$gr(pb$opt$par), digits = 4), "\n")
cat("  warnings body:", attr(pb, "wc"), " lin:", attr(pl, "wc"), "\n")
# the same fit with the body coefficient started at the linear optimum
pb2 <- tt(frm(bf(y ~ a * xt, a ~ 1, nl = TRUE), data = d,
              prior = prior(normal(0, 1e-3), nlpar = a),
              start = list(beta = c("a_(Intercept)" = 0),
                           betad = pl$estimates$betad)))
cat(sprintf("  body started at lin's sigma: logLik %.8f sigma %.10f\n",
            llv(pb2), pb2$estimates$betad[[1]]))
cat("=== prior normal(5e19, 1e19) on a (located), tiny x\n")
pb <- tt(frm(bf(y ~ a * xt, a ~ 1, nl = TRUE), data = d,
             prior = prior(normal(5e19, 1e19), nlpar = a)))
pl <- tt(frm(y ~ 0 + xt, data = d,
             prior = prior(normal(5e19, 1e19), class = b)))
s <- exp(pb$estimates$betad[[1]])
amap <- (sum(d$xt * d$y) / s^2 + 5e19 / 1e19^2) /
  (sum(d$xt^2) / s^2 + 1 / 1e19^2)
cat(sprintf("body a %.8g  lin b %.8g  exact MAP %.8g  rel body %.2e lin %.2e\n",
            pb$estimates$beta[[1]], pl$estimates$beta[[1]], amap,
            abs(pb$estimates$beta[[1]] / amap - 1),
            abs(pl$estimates$beta[[1]] / amap - 1)))
cat("  se body", sqrt(diag(vcov(pb)))[1], "lin", sqrt(diag(vcov(pl)))[1], "\n")

cat("=== user start for a scaled coefficient\n")
for (a0 in c(0, 1, 5e19, -1e21)) {
  f <- tt(frm(bf(y ~ a * xt, a ~ 1, nl = TRUE), data = d,
              start = list(beta = c("a_(Intercept)" = a0))))
  cat(sprintf("start a = %g: logLik %.8f (lin %.8f) a %.6g conv %s units %s",
              a0, llv(f), llv(fl), f$estimates$beta[[1]], f$opt$convergence,
              paste(format(f$par_units %||% 1, digits = 4), collapse = ",")),
      "start echo:", format(unlist(f$start)), "\n")
}

cat("=== large y (seed 9901): a + b * x\n")
set.seed(9901)
d2 <- data.frame(x = rnorm(n))
d2$y <- 1e6 + 5e4 * d2$x + rnorm(n, 0, 1e3)
fb2 <- tt(frm(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), data = d2))
fl2 <- tt(frm(y ~ x, data = d2))
show("a + b * x vs y ~ x", fb2, fl2, c("a_(Intercept)", "b_(Intercept)"),
     c("(Intercept)", "x"))

cat("=== added-term detection and units per body\n")
set.seed(4242)
d3 <- data.frame(x = rnorm(n), xp = runif(n, 0.5, 3))
d3$y <- 1 + 0.5 * d3$x + rnorm(n)
bodies <- list(
  "a + b * x" = bf(y ~ a + b * x, a + b ~ 1, nl = TRUE),
  "(a + b) * x" = bf(y ~ (a + b) * x, a + b ~ 1, nl = TRUE),
  "a * x + c" = bf(y ~ a * x + c, a + c ~ 1, nl = TRUE),
  "exp(a) + b" = bf(y ~ exp(a) + b, a + b ~ 1, nl = TRUE),
  "a - b * x" = bf(y ~ a - b * x, a + b ~ 1, nl = TRUE),
  "a * b * x" = bf(y ~ a * b * x, a + b ~ 1, nl = TRUE),
  "a / xp + b" = bf(y ~ a / xp + b, a + b ~ 1, nl = TRUE),
  "a * x^2 + b * x" = bf(y ~ a * x^2 + b * x, a + b ~ 1, nl = TRUE),
  "a * (x > 0) + b" = bf(y ~ a * (x > 0) + b, a + b ~ 1, nl = TRUE),
  "b1 * exp(b2 * x)" = bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
  "a^2 * x + b" = bf(y ~ a^2 * x + b, a + b ~ 1, nl = TRUE),
  "a * x, a ~ 1 + xp" = bf(y ~ a * x, a ~ 1 + xp, nl = TRUE))
for (nm in names(bodies)) {
  f <- tt(withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
          frm(bodies[[nm]], data = d3)))
  if (iserr(f)) { cat(nm, ": ERR", f$msg, "\n"); next }
  fr <- f$frame
  lp <- Filter(function(l) !is.null(l[["nl_body"]]), fr[["linpreds"]])[[1]]
  bc <- frmtmb:::nl_body_coefs(fr, lp)$coefs
  add <- frmtmb:::nl_body_additive(fr, lp, f$estimates, bc)
  cat(sprintf("%-20s coefs %-28s added %s\n", nm,
              paste(names(unlist(bc)), collapse = ","),
              paste(add, collapse = ",")))
}
