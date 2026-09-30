# Re-check (punch round 1): the new refusals with their guard-absent
# cases (item 3), formula-environment constants (item 4) and the stored
# call of a pooled update() (item 5). Data seed 93. REVLIB="" for base.
LIB <- Sys.getenv("REVLIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
tr <- function(label, expr) {
  w <- character(0); m <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }, message = function(c) {
    m <<- c(m, trimws(conditionMessage(c))); invokeRestart("muffleMessage")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  cat(sprintf("[%s] %s\n", label, if (inherits(r, "err"))
    paste("ERROR:", substr(r, 1, 260)) else "ok"))
  if (length(w)) cat(sprintf("   WARN: %s\n", substr(unique(w), 1, 200)))
  if (length(m)) cat(sprintf("   MSG: %s\n", substr(unique(m), 1, 200)))
  invisible(if (inherits(r, "err")) NULL else r)
}
same <- function(label, a, b) cat(sprintf("   [%s] identical: %s\n", label,
                                          identical(a, b)))
bf <- frmtmb::bf
set.seed(93)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n), w = runif(n, 0.5, 2),
                nt = rpois(n, 5) + 3L, s = runif(n, 0.2, 0.5),
                t = rep(1:6, 10), g = factor(rep(1:10, each = 6)))
d$y <- rnorm(n, 0.5 * d$x)
d$yb <- rbinom(n, d$nt, 0.4)
d$y01 <- rbinom(n, 1, 0.5)

cat("\n== 3a. NA refusal: guard-absent cases (NA in response, predictor, expression variable) ==\n")
for (arm in c("response NA", "predictor NA", "expr variable NA",
              "expr NA", "mi() + expr variable NA")) {
  dd <- d
  f <- bf(y | weights(w * 2) ~ x)
  if (arm == "response NA") dd$y[c(3, 9)] <- NA
  if (arm == "predictor NA") dd$x[c(3, 9)] <- NA
  if (arm == "expr variable NA") dd$w[c(3, 9)] <- NA
  if (arm == "expr NA") f <- bf(y | weights(ifelse(x > 1.5, NA, w)) ~ x)
  if (arm == "mi() + expr variable NA") {
    dd$xm <- dd$x; dd$xm[c(4, 10)] <- NA; dd$w[c(3, 9)] <- NA
    f <- bf(y | weights(w * 2) ~ mi(xm)) + bf(xm | mi() ~ z)
  }
  r <- tr(paste("weights(w * 2),", arm), frm(f, data = dd))
  if (!is.null(r)) cat("   nobs", nobs(r), "\n")
  if (arm != "expr NA" && arm != "mi() + expr variable NA") {
    dd$w2 <- dd$w * 2
    r2 <- tr(paste("weights(w2),", arm), frm(bf(y | weights(w2) ~ x),
                                            data = dd))
    if (!is.null(r) && !is.null(r2)) same("expr vs column logLik",
                                          logLik(r), logLik(r2))
  }
}
d$sub <- d$x > 0
dd <- d; dd$w[!dd$sub][1:3] <- NA
r <- tr("subset(): expression NA only on rows outside the subset",
        frm(bf(y | weights(ifelse(is.na(w), NA, w * 2)) + subset(sub) ~ x) +
              bf(yb | trials(nt) ~ x, family = binomial()), data = dd))
dd$ex <- ifelse(dd$x > 0, dd$w, NA)
r <- tr("subset(): expression NA from complete vars only outside subset",
        frm(bf(y | weights(ifelse(x > 0, w, NA)) + subset(sub) ~ x) +
              bf(yb | trials(nt) ~ x, family = binomial()), data = d))
tr("brms same", suppressWarnings(brms::standata(
  brms::bf(y | weights(ifelse(x > 0, w, NA)) + subset(sub) ~ x) +
    brms::bf(yb | trials(nt) ~ x, family = stats::binomial()) +
    brms::set_rescor(FALSE), data = d)))
dc <- d
dc$cc <- ifelse(seq_len(n) %% 3 == 0, "interval", "none")
dc$hi <- ifelse(dc$cc == "interval", dc$y + 1, NA)
tr("cens interval bound NA on non-interval rows (expression)",
   frm(bf(y | cens(cc, hi * 1) ~ x), data = dc))
tr("mi() response NA with weights(w * 2)", {
  dm <- d; dm$y[c(2, 5)] <- NA
  frm(bf(y | mi() + weights(w * 2) ~ x), data = dm) })
tr("vapply(w, round, 1): a function as an argument",
   frm(bf(y | weights(vapply(w, round, 1)) ~ x), data = d))
tr("brms vapply(w, round, 1)", suppressWarnings(brms::standata(
  y | weights(vapply(w, round, 1)) ~ x, data = d)))

cat("\n== 3b. refit() refusal ==\n")
fb <- frm(bf(y01 ~ x), data = d, family = bernoulli())
tr("refit 0/1", refit(fb, simulate(fb, seed = 1)[[1]]))
tr("refit logical", refit(fb, simulate(fb, seed = 1)[[1]] == 1))
tr("refit with NA", { r <- simulate(fb, seed = 1)[[1]]; r[2] <- NA
  refit(fb, r) })
tr("refit value 2", refit(fb, rep(c(0, 2), n / 2)))
fg <- frm(bf(y ~ x), data = d)
tr("gaussian refit any values", refit(fg, d$y * 10))
fbin <- frm(bf(yb | trials(nt) ~ x), data = d, family = binomial())
tr("binomial refit counts", refit(fbin, pmin(d$yb + 1L, d$nt)))
d$ych <- ifelse(d$y01 == 1, "yes", "no")
fc <- frm(bf(ych ~ x), data = d, family = bernoulli())
tr("character-response refit with 0/1", refit(fc, d$y01))
tr("character-response refit with 'yes'/'no'", refit(fc, d$ych))
tr("0/1 fit, mixture(bernoulli)?", NULL)

cat("\n== 3c. single-bound cens(), weights(t * 2), mv update ==\n")
dc$ub <- max(dc$y) + 10
tr("cens(cc, ub) column", frm(bf(y | cens(cc, ub) ~ x), data = dc))
tr("cens(cc, 10) literal", frm(bf(y | cens(cc, 10) ~ x), data = dc))
tr("cens(cc, max(y) + 10)", frm(bf(y | cens(cc, max(y) + 10) ~ x),
                                 data = dc))
tr("cens(cc) no bound, no interval rows", {
  d2 <- dc; d2$cc <- "none"; frm(bf(y | cens(cc) ~ x), data = d2) })
tr("weights(t * 2) with column t", frm(bf(y | weights(t * 2) ~ x),
                                       data = d))
d_not <- d; d_not$t <- NULL
tr("weights(t * 2) without column t", frm(bf(y | weights(t * 2) ~ x),
                                          data = d_not))
tvec <- rep(1:6, 10)
f_tv <- local({ t <- tvec; y | weights(t * 2) ~ x })
tr("weights(t * 2), t a vector of the formula env", frm(bf(f_tv),
                                                        data = d_not))
tr("offset(log(t)) without column t", frm(bf(y ~ x + offset(log(t))),
                                          data = d_not))
fm <- frm(bf(y ~ x, sigma ~ z) + bf(yb | trials(nt) ~ x,
                                    family = binomial()), data = d)
tr("mv update complete formula", update(fm, bf(y ~ z) +
                                          bf(yb | trials(nt) ~ x,
                                             family = binomial())))
tr("mv update newdata only", update(fm, newdata = d[1:50, ]))
tr("mv update delta", update(fm, ~ . + z))
tr("mv update, no formula, iter control", update(fm, control =
                                                   frmtmb_control()))
fu <- frm(bf(y ~ x, sigma ~ z), data = d)
tr("univariate -> mv update", update(fu, bf(y ~ x) +
                                       bf(yb | trials(nt) ~ x,
                                          family = binomial())))

cat("\n== 4. formula-environment constants ==\n")
k <- 12L
d$yk <- rbinom(n, 12, 0.3)
ck <- 0.5
f_tk <- frm(bf(yk | trials(k) ~ x), data = d, family = binomial())
f_tl <- frm(bf(yk | trials(12) ~ x), data = d, family = binomial())
same("trials(k) vs trials(12) logLik", logLik(f_tk), logLik(f_tl))
f_pi <- tr("plain predictor I(x * ck)", frm(bf(y ~ I(x * ck)), data = d))
f_off <- tr("offset(z * ck)", frm(bf(y ~ x + offset(z * ck)), data = d))
f_dp <- tr("dpar formula sigma ~ I(z * ck)", frm(bf(y ~ x, sigma ~ I(z * ck)),
                                                 data = d))
f_w <- frm(bf(y | weights(w * ck) ~ x), data = d)
nd <- d[1:4, ]
before <- list(tk = fitted(f_tk, newdata = nd)[, 1],
               pi = fitted(f_pi, newdata = nd)[, 1],
               off = fitted(f_off, newdata = nd)[, 1],
               dp = fitted(f_dp, newdata = nd, dpar = "sigma")[, 1],
               ll_tk = logLik(f_tk), ll_pi = logLik(f_pi))
k <- 20L; ck <- 2
after <- list(tk = fitted(f_tk, newdata = nd)[, 1],
              pi = fitted(f_pi, newdata = nd)[, 1],
              off = fitted(f_off, newdata = nd)[, 1],
              dp = fitted(f_dp, newdata = nd, dpar = "sigma")[, 1])
for (nm in c("tk", "pi", "off", "dp")) {
  cat(sprintf("   after k/ck change, fitted(newdata) %s unchanged: %s\n",
              nm, identical(before[[nm]], after[[nm]])))
}
cat("   in-sample fitted of trials(k) still the fit's:",
    identical(fitted(f_tk)[, 1], fitted(f_tl)[, 1]), "\n")
u_tk <- tr("update(f_tk) after k change", update(f_tk))
if (!is.null(u_tk)) cat("   update logLik unchanged:",
                        identical(logLik(u_tk), before$ll_tk), "\n")
u_pi <- tr("update(f_pi) after ck change", update(f_pi))
if (!is.null(u_pi)) cat("   update I(x*ck) coef ratio:",
                        coef(u_pi)[2] / coef(f_pi)[2], "\n")
tr("predict(f_tk, newdata) after k change", {
  set.seed(1); p <- predict(f_tk, newdata = nd)
  cat("   predict max after k=20:", max(p[, "Q97.5"]), "\n") })
tr("simulate(f_tk) in sample after k change", {
  s <- simulate(f_tk, nsim = 20, seed = 1)
  cat("   simulate max:", max(unlist(s)), "(k was 12 at fit)\n") })
tr("residuals(f_tk) in sample after k change",
   same("residuals trials(k) vs trials(12)", residuals(f_tk),
        residuals(f_tl)))
tr("influence(f_tk) after k change", {
  a <- influence(frm(bf(yk | trials(12) ~ x + (1 | g)), data = d,
                     family = binomial()), groups = "g")$fixed
  k <- 12L
  f2 <- frm(bf(yk | trials(k) ~ x + (1 | g)), data = d, family = binomial())
  k <<- 20L
  b <- influence(f2, groups = "g")$fixed
  cat("   influence refits of trials(k) after k change equal trials(12):",
      identical(a, b), "\n") })

cat("\n== 5. eval(fit$call) after a pooled update() ==\n")
k <- 12L; ck <- 0.5
d$yo <- cut(d$y, c(-Inf, -0.5, 0, 0.5, Inf), labels = FALSE)
d$ycat <- factor(sample(c("p", "q", "r"), n, TRUE))
d$ypos <- exp(rnorm(n, 0.2 * d$x, 0.3))
cases <- list(
  dpar = function() frm(bf(y ~ x, sigma ~ z), data = d),
  student_link = function() frm(bf(y ~ x, sigma ~ z,
                                   family = student(link_sigma = "softplus")),
                                data = d),
  cumul_equid = function() frm(bf(yo ~ x, disc ~ z,
                                  family = cumulative(threshold =
                                                        "equidistant")),
                               data = d),
  cumul_equid_plus = function() frm(bf(yo ~ x, disc ~ z) +
                                      cumulative(threshold = "equidistant"),
                                    data = d),
  cumul_equid_arg = function() frm(bf(yo ~ x, disc ~ z), data = d,
                                   family = cumulative(threshold =
                                                         "equidistant")),
  categ_refcat = function() frm(bf(ycat ~ x, muq ~ z,
                                   family = categorical(refcat = "r")),
                                data = d),
  gamma_link = function() frm(bf(ypos ~ x, shape ~ z,
                                 family = Gamma(link = "identity")),
                              data = d),
  sigma_const = function() frm(bf(y ~ x, sigma = 1.5), data = d),
  nl = function() frm(bf(y ~ a + b * x, a ~ 1 + z, b ~ 1, nl = TRUE),
                      data = d))
newf <- list(dpar = y ~ x + z, student_link = y ~ x + z,
             cumul_equid = yo ~ x + z, cumul_equid_plus = yo ~ x + z,
             cumul_equid_arg = yo ~ x + z, categ_refcat = ycat ~ x + z,
             gamma_link = ypos ~ x + z, sigma_const = y ~ x + z,
             nl = y ~ a + b * x)
for (nm in names(cases)) {
  f0 <- tr(paste(nm, "fit"), cases[[nm]]())
  if (is.null(f0)) next
  u <- tr(paste(nm, "update"), update(f0, newf[[nm]]))
  if (is.null(u)) next
  e <- tr(paste(nm, "eval(u$call)"), eval(u$call))
  if (is.null(e)) next
  same(paste(nm, "eval(call) logLik"), logLik(e), logLik(u))
  same(paste(nm, "eval(call) coef names"), rownames(fixef(e)),
       rownames(fixef(u)))
  fam0 <- f0$spec$responses[[1]]$family; famu <- u$spec$responses[[1]]$family
  for (fld in c("threshold", "refcat")) {
    if (!is.null(fam0[[fld]]) || !is.null(famu[[fld]])) {
      cat(sprintf("   %s: original %s, updated %s\n", fld,
                  format(fam0[[fld]]), format(famu[[fld]])))
    }
  }
  cat(sprintf("   n pars original %d, updated %d; call: %s\n",
              length(f0$obj$par), length(u$obj$par),
              substr(deparse1(u$call$formula), 1, 160)))
}
