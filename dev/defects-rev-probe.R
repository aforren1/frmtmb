# Reviewer of lane defects: behavioral probes of each fix on the lane
# build, with the absent case of each guard. Seeds are set per block.
#   Rscript dev/defects-rev-probe.R
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
hdr <- function(...) cat("\n=====", ..., "\n")
try_msg <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  if (length(w)) cat("  [warnings:", paste(unique(w), collapse = " | "), "]\n")
  v
}
say <- function(label, v) {
  if (is.character(v) && length(v) == 1) cat(" ", label, ":", substr(v, 1, 300), "\n")
  else { cat(" ", label, ":\n"); print(v) }
}

## 1. cs() linear layers: brms convention and the standard error --------
hdr("1. fitted(scale = 'linear') on cs() fits")
set.seed(11)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, 20)))
lat <- 0.6 * d$x + rnorm(n)
d$y <- cut(lat, c(-Inf, -0.7, 0.2, 1, Inf), labels = FALSE)
fit <- frm(y ~ x + cs(z), family = sratio(), data = d)
L <- fitted(fit, scale = "linear")
say("dim", dim(L)); say("dimnames", dimnames(L))
est <- fit$estimates
say("estimate names", names(est))
X <- model.matrix(~ x, d)
# the shared predictor (no intercept: thresholds carry it) and the cs part
b <- fixef(fit)
say("fixef rows", rownames(b))
# independent delta-method SE of layer k at row i, from the joint
# covariance of the outer parameters via numerical Jacobian
lin_fun <- function(par) {
  fit2 <- fit
  fit2$opt$par <- par
  NULL
}
eta_k <- L[, "Estimate", ]
# identity check: sratio P(Y = 1) = plogis(tau_1 - eta_1)
P <- fitted(fit)
tau <- fit$estimates[[grep("thres|Intercept", names(fit$estimates))[1]]]
say("thresholds", tau)
p1 <- plogis(tau[1] - eta_k[, 1])
say("max |P(Y=1) - F(tau1 - eta1)| / max P", max(abs(P[, "Estimate", 1] - p1)) / max(P[, "Estimate", 1]))
# brms convention: .predictor_cs adds b[, I + k] %*% X to eta, and the
# sratio density reads thres - eta; so layer k minus layer j equals
# (bcs_k - bcs_j) * z
cs_est <- est[[grep("cs", names(est))[1]]]
say("cs coefficient vector", cs_est)
dd <- (eta_k[, 2] - eta_k[, 1]) / d$z
say("range of (eta_2 - eta_1)/z", range(dd))
say("cs_2 - cs_1", cs_est[2] - cs_est[1])
# SE: compare with a hand delta method through fit$sdr covariance
se_fd <- L[1:3, "Est.Error", ]
say("Est.Error rows 1:3", se_fd)
se_shared <- try_msg(fitted(frm(y ~ x + z, family = sratio(), data = d),
                            scale = "linear")[1:3, "Est.Error"])
say("Est.Error without cs (shared eta only)", se_shared)
# newdata path and one row
nd <- d[1:3, ]
Ln <- try_msg(fitted(fit, newdata = nd, scale = "linear"))
say("newdata layers equal in-sample rows", isTRUE(all.equal(unclass(Ln[, "Estimate", ]), unclass(L[1:3, "Estimate", ]), check.attributes = FALSE)))
Ln1 <- try_msg(fitted(fit, newdata = d[1, , drop = FALSE], scale = "linear"))
say("one-row newdata dim", dim(Ln1))
# response scale unchanged, dpar = "mu"
say("fitted(dpar = 'mu', scale = 'linear') dim", dim(try_msg(fitted(fit, dpar = "mu", scale = "linear"))))
say("frm_linpred link dim (unchanged API)", dim(try_msg(frm_linpred(fit, type = "link"))))
# na.exclude: in-sample padding
d_na <- d; d_na$x[5] <- NA
fit_na <- frm(y ~ x + cs(z), family = sratio(), data = d_na, na.action = na.exclude)
Lna <- try_msg(fitted(fit_na, scale = "linear"))
say("na.exclude dim", dim(Lna))
say("row 5 NA in every layer", all(is.na(Lna[5, "Estimate", ])))
# cs with random effects
fit_re <- frm(y ~ cs(x) + (1 | g), family = acat(), data = d)
Lre <- try_msg(fitted(fit_re, scale = "linear"))
say("acat cs + re dim", dim(Lre))
Lre0 <- try_msg(fitted(fit_re, scale = "linear", re_formula = NA))
say("re_formula = NA dim", dim(Lre0))

## 2. residuals: newdata, multivariate --------------------------------
hdr("2. residuals(newdata =), multivariate residuals")
set.seed(21)
n <- 100
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, 10)),
                tr = 8L)
d$y <- 1 + d$x + rnorm(10)[d$g] + rnorm(n)
d$pos <- exp(0.2 + 0.3 * d$x + rnorm(n, 0, 0.4))
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x))
d$bin <- rbinom(n, d$tr, plogis(d$x))
d$y2 <- 0.5 * d$y + rnorm(n)
chk <- function(label, fit) {
  r0 <- try_msg(residuals(fit))
  r1 <- try_msg(residuals(fit, newdata = fit$data))
  if (is.character(r0) || is.character(r1)) { say(label, if (is.character(r1)) r1 else r0); return() }
  cat(sprintf("  %-12s in-sample vs newdata=data: max|dEst| = %.3g, max|dSE| = %.3g\n",
              label, max(abs(r0[, 1] - r1[, 1])), max(abs(r0[, 2] - r1[, 2]))))
}
chk("gaussian", frm(y ~ x, data = d))
chk("gauss re", frm(y ~ x + (1 | g), data = d))
chk("lognormal", frm(pos ~ x, family = lognormal(), data = d))
chk("log(y)", frm(log(pos) ~ x, data = d))
chk("poisson", frm(cnt ~ x, family = poisson(), data = d))
chk("binomial", frm(bin | trials(tr) ~ x, family = binomial(), data = d))
chk("dist", frm(bf(y ~ x, sigma ~ z), data = d))
# new levels
fre <- frm(y ~ x + (1 | g), data = d)
ndn <- data.frame(x = c(0, 1), g = factor(c("new1", "1")), y = c(1, 2))
say("newdata new level, allow_new_levels = FALSE", try_msg(residuals(fre, newdata = ndn)))
rn <- try_msg(residuals(fre, newdata = ndn, allow_new_levels = TRUE))
say("newdata new level, allow = TRUE", rn)
fn <- try_msg(fitted(fre, newdata = ndn, allow_new_levels = TRUE))
say("equals y - fitted", isTRUE(all.equal(unname(rn[, 1]), ndn$y - unname(fn[, 1]))))
say("pearson + newdata", try_msg(residuals(fre, type = "pearson", newdata = ndn)))
say("newdata without response", try_msg(residuals(fre, newdata = ndn[, 1:2])))
say("positional newdata (2nd arg is type)", try_msg(residuals(fre, ndn)))
# multivariate
fmv <- frm(bf(y ~ x) + bf(y2 ~ z) + set_rescor(TRUE), data = d)
R <- try_msg(residuals(fmv))
say("mv dim", dim(R)); say("mv dimnames[[3]]", dimnames(R)[[3]])
R1 <- try_msg(residuals(fmv, resp = "y2"))
say("resp = 'y2' dim", dim(R1))
say("layer y2 equals resp = 'y2'", isTRUE(all.equal(unclass(R[, , "y2"]), unclass(R1), check.attributes = FALSE)))
F2 <- fitted(fmv, resp = "y2")
say("layer y2 equals y2 - fitted", isTRUE(all.equal(unname(R[, "Estimate", "y2"]), d$y2 - unname(F2[, "Estimate"]))))
say("resp = c('y2','y') order", dimnames(try_msg(residuals(fmv, resp = c("y2", "y"))))[[3]])
say("resp unknown", try_msg(residuals(fmv, resp = "nope")))
say("mv osa", try_msg(residuals(fmv, type = "osa")))
say("mv newdata", dim(try_msg(residuals(fmv, newdata = d[1:5, ]))))
# mv with mi(): NA rows of the imputed response
d$xmi <- ifelse(runif(n) < 0.2, NA, d$x)
fmi <- frm(bf(y ~ mi(xmi)) + bf(xmi | mi() ~ z) + set_rescor(FALSE), data = d)
Rmi <- try_msg(residuals(fmi))
say("mi dim", dim(Rmi))
say("xmi layer NA where missing", if (is.array(Rmi)) table(is.na(Rmi[, "Estimate", "xmi"]), is.na(d$xmi)) else Rmi)
# mv with a response dropped on NA rows (non-mi): are rows aligned?
d$y3 <- d$y2; d$y3[c(2, 7)] <- NA
fna <- try_msg(frm(bf(y ~ x) + bf(y3 ~ z) + set_rescor(FALSE), data = d))
if (!is.character(fna)) {
  Rna <- try_msg(residuals(fna))
  say("NA-row mv dim", dim(Rna))
  say("nobs", nobs(fna))
}
fna2 <- try_msg(frm(bf(y ~ x) + bf(y3 ~ z) + set_rescor(FALSE), data = d,
                    na.action = na.exclude))
if (!is.character(fna2)) {
  Rna2 <- try_msg(residuals(fna2))
  say("na.exclude mv dim", dim(Rna2))
  say("per resp dims", sapply(c("y", "y3"), function(r) nrow(residuals(fna2, resp = r))))
}

## 3. predict(sample_new_levels = "old_levels") ------------------------
hdr("3. predict(sample_new_levels = 'old_levels')")
set.seed(31)
n <- 200
d <- data.frame(g = factor(rep(1:10, 20)), x = rnorm(n))
u <- c(-6, -4, -2, 0, 1, 2, 3, 4, 5, 6)
d$y <- u[d$g] + rnorm(n, 0, 0.3)
f <- frm(y ~ 1 + (1 | g), data = d)
re <- ranef(f)$g[, 1] + fixef(f)[1, 1]
nd <- data.frame(g = factor(c("A", "A", "B")), x = 0)
set.seed(1)
p <- predict(f, newdata = nd, allow_new_levels = TRUE,
             sample_new_levels = "old_levels", ndraws = 2000, summary = FALSE)
say("draws dim", dim(p))
m <- colMeans(p)
say("column means", m)
say("nearest seen level per column", sapply(m, function(v) which.min(abs(re - v))))
say("rows of same new level share a level (A,A)", isTRUE(all.equal(m[1], m[2], tolerance = 0.05)))
say("sd per column (sigma ~ 0.3)", apply(p, 2, sd))
set.seed(1)
p2 <- predict(f, newdata = nd, allow_new_levels = TRUE,
              sample_new_levels = "old_levels", ndraws = 2000, summary = FALSE)
say("reproducible under set.seed", identical(p, p2))
# how often do A and B pick different levels over calls
picks <- replicate(40, {
  pp <- predict(f, newdata = nd, allow_new_levels = TRUE,
                sample_new_levels = "old_levels", ndraws = 20, summary = FALSE)
  sapply(colMeans(pp), function(v) which.min(abs(re - v)))
})
say("distinct picks for A over 40 calls", length(unique(picks[1, ])))
say("A != B in calls", mean(picks[1, ] != picks[3, ]))
say("gaussian sd per column", apply(predict(f, newdata = nd, allow_new_levels = TRUE, ndraws = 2000, summary = FALSE), 2, sd))
say("uncertainty refused", try_msg(predict(f, newdata = nd, allow_new_levels = TRUE, sample_new_levels = "uncertainty")))
say("fitted old_levels refused", try_msg(fitted(f, newdata = nd, allow_new_levels = TRUE, sample_new_levels = "old_levels")))
say("fitted gaussian == default", identical(fitted(f, newdata = nd, allow_new_levels = TRUE, sample_new_levels = "gaussian"), fitted(f, newdata = nd, allow_new_levels = TRUE)))
say("sample_new_levels bad value", try_msg(predict(f, newdata = nd, allow_new_levels = TRUE, sample_new_levels = "gaus")))
# correlated slopes: the picked level's intercept AND slope together
set.seed(32)
d$y <- u[d$g] + (u[d$g] / 3) * d$x + rnorm(n, 0, 0.3)
f2 <- frm(y ~ x + (1 + x | g), data = d)
nd2 <- data.frame(g = factor(c("A", "A")), x = c(0, 1))
set.seed(5)
pp <- predict(f2, newdata = nd2, allow_new_levels = TRUE,
              sample_new_levels = "old_levels", ndraws = 3000, summary = FALSE)
co <- coef(f2)$g
lev <- which.min(abs(co[, "Estimate", "Intercept"] - mean(pp[, 1])))
say("picked level", lev)
say("slope implied vs that level's slope", c(mean(pp[, 2]) - mean(pp[, 1]), co[lev, "Estimate", "x"]))
# sigma dpar random effect and old_levels
set.seed(33)
d$y <- u[d$g] + rnorm(n, 0, exp(0.3 * (as.integer(d$g) - 5) / 5))
f3 <- try_msg(frm(bf(y ~ 1 + (1 | g), sigma ~ 1 + (1 | g)), data = d))
if (!is.character(f3)) {
  say("dpar re old_levels", dim(try_msg(predict(f3, newdata = nd, allow_new_levels = TRUE, sample_new_levels = "old_levels", ndraws = 50, summary = FALSE))))
}
# by = groups: brms picks from levels with the same `by`
d$grp <- factor(ifelse(as.integer(d$g) <= 5, "lo", "hi"))
f4 <- try_msg(frm(y ~ 1 + (1 | gr(g, by = grp)), data = d))
if (!is.character(f4)) {
  nd4 <- data.frame(g = factor("A"), grp = factor("lo", levels = c("hi", "lo")))
  set.seed(9)
  pp4 <- try_msg(replicate(20, colMeans(predict(f4, newdata = nd4, allow_new_levels = TRUE, sample_new_levels = "old_levels", ndraws = 30, summary = FALSE))))
  say("by = grp, new level in 'lo': picked levels (lo are 1..5, means -6..0)", if (is.character(pp4)) pp4 else sapply(pp4, function(v) which.min(abs(u - v))))
}

## 4. newdata factor coded as numbers, contrasts attribute -------------
hdr("4. numeric codes for a factor, contrasts attribute")
set.seed(41)
n <- 80
d <- data.frame(Trt = factor(sample(c(0, 1), n, TRUE)), x = rnorm(n))
d$y <- 1 + as.numeric(as.character(d$Trt)) + d$x + rnorm(n)
f <- frm(y ~ Trt + x, data = d)
a1 <- fitted(f, newdata = data.frame(Trt = c(0, 1), x = 0))
a2 <- fitted(f, newdata = data.frame(Trt = factor(c("0", "1")), x = 0))
say("numeric codes == factor", identical(a1, a2))
say("code not a level", try_msg(fitted(f, newdata = data.frame(Trt = 2, x = 0))))
say("code 1.0 (double)", identical(fitted(f, newdata = data.frame(Trt = 1.0, x = 0)), a2[2, , drop = FALSE]))
# character labels of a factor with levels "a","b"
d$f2 <- factor(sample(c("a", "b"), n, TRUE))
f5 <- frm(y ~ f2 + x, data = d)
say("numeric in a letter-level factor", try_msg(fitted(f5, newdata = data.frame(f2 = 1, x = 0))))
# a factor fitted with levels 1:3 and newdata the integer 2
d$k <- factor(sample(1:3, n, TRUE))
f6 <- frm(y ~ k + x, data = d)
say("integer code 2", identical(fitted(f6, newdata = data.frame(k = 2L, x = 0)), fitted(f6, newdata = data.frame(k = factor("2", levels = 1:3), x = 0))))
# logical predictor fitted as logical
d$lg <- d$x > 0
f7 <- frm(y ~ lg, data = d)
say("logical predictor newdata", try_msg(fitted(f7, newdata = data.frame(lg = c(TRUE, FALSE)))))
# a factor with a contrasts attribute: slice predictions, no warning
d2 <- d; contrasts(d2$Trt) <- contr.sum(2)
f8 <- frm(y ~ Trt + x, data = d2)
w8 <- try_msg(fitted(f8, newdata = d2[1:10, ]))
say("contr.sum fit: slice equals in-sample", isTRUE(all.equal(unclass(w8), unclass(fitted(f8)[1:10, ]), check.attributes = FALSE)))
# newdata carrying a DIFFERENT contrasts attribute: the fitted coding wins
d3 <- d2[1:10, ]; contrasts(d3$Trt) <- contr.treatment(2)
w9 <- try_msg(fitted(f8, newdata = d3))
say("newdata with other contrasts equals fitted coding", isTRUE(all.equal(unclass(w9), unclass(w8), check.attributes = FALSE)))
# a factor in a group-level term coded as a number
d$gg <- factor(rep(1:8, 10))
f9 <- frm(y ~ x + (1 | gg), data = d)
say("numeric group code", identical(fitted(f9, newdata = data.frame(x = 0, gg = 3)), fitted(f9, newdata = data.frame(x = 0, gg = factor("3", levels = 1:8)))))
say("numeric new group code, not allowed", try_msg(fitted(f9, newdata = data.frame(x = 0, gg = 99))))

## 5. formula helpers inside a formula ---------------------------------
hdr("5. helpers inside a formula, absent cases")
set.seed(51)
d <- data.frame(y = rnorm(40), x = rnorm(40), y2 = rnorm(40))
say("set_rescor in rhs", try_msg(frm(y ~ x + set_rescor(TRUE), data = d)))
say("lf in rhs", try_msg(frm(y ~ x + lf(sigma ~ x), data = d)))
lf <- function(v) v^2
say("user's own lf() is a term", try_msg(fixef(frm(y ~ lf(x), data = d))))
rm(lf)
lf_col <- try_msg(fixef(frm(y ~ x, data = transform(d, lf = x))))
say("a column named lf", lf_col)
# unattached frmtmb: frmtmb::frm with no helper visible from globalenv
res_unatt <- local({
  e <- new.env(parent = baseenv())
  e$d <- d
  eval(quote(tryCatch(fixef(frmtmb::frm(frmtmb::bf(y ~ x), data = d)),
                      error = function(err) conditionMessage(err))), e)
})
say("frmtmb::frm with nothing attached (no helper in formula)", res_unatt)
res_unatt2 <- local({
  e <- new.env(parent = baseenv()); e$d <- d
  eval(quote(tryCatch(frmtmb::frm(y ~ x + set_rescor(TRUE), data = d),
                      error = function(err) conditionMessage(err))), e)
})
say("frmtmb::frm unattached with set_rescor in rhs", res_unatt2)
say("nlf body y ~ nlf(...)? nl formula legit", try_msg(fixef(frm(bf(y ~ a * x, a ~ 1, nl = TRUE), data = d))))
say("mv set_rescor legit", try_msg(class(frm(bf(y ~ x) + bf(y2 ~ x) + set_rescor(FALSE), data = d))))
say("lf legit", try_msg(class(frm(bf(y ~ x) + lf(sigma ~ x), data = d))))

## 6. cs() in a bar, absent cases -------------------------------------
hdr("6. (cs(x) | g)")
set.seed(61)
n <- 120
d <- data.frame(x = rnorm(n), g = factor(rep(1:12, 10)))
d$y <- cut(d$x + rnorm(n), c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
d$yc <- rnorm(n)
say("sratio (cs(x) | g)", try_msg(frm(y ~ x + (cs(x) | g), family = sratio(), data = d)))
say("cumulative (cs(x) | g)", try_msg(frm(y ~ x + (cs(x) | g), family = cumulative(), data = d)))
say("gaussian (cs(x) | g)", try_msg(frm(yc ~ x + (cs(x) | g), data = d)))
say("sratio cs(x) + (1 | g) legit", try_msg(class(frm(y ~ cs(x) + (1 | g), family = sratio(), data = d))))
say("cs(1 | g) compound symmetry legit", try_msg(class(frm(yc ~ x + cs(1 | g), data = d))))
say("cs(1 + x | g) compound symmetry legit", try_msg(class(frm(yc ~ x + cs(1 + x | g), data = d))))
say("sratio + cs(1 | g) covariance", try_msg(class(frm(y ~ x + cs(1 | g), family = sratio(), data = d))))

## 7. se() early check, absent cases ----------------------------------
hdr("7. se()")
set.seed(71)
d <- data.frame(y = rnorm(50), x = rnorm(50), s = runif(50, 0.2, 1))
say("gaussian se present", try_msg(class(frm(y | se(s) ~ x, data = d))))
say("gaussian se column missing", try_msg(frm(y | se(nope) ~ x, data = d)))
say("weibull se column missing", try_msg(frm(abs(y) | se(nope) ~ x, family = weibull(), data = d)))
say("weibull se present", try_msg(frm(abs(y) | se(s) ~ x, family = weibull(), data = d)))
sfree <- runif(50, 0.2, 1)
say("gaussian se from formula env", try_msg(class(frm(y | se(sfree) ~ x, data = d))))
say("poisson se from env missing col", try_msg(frm(rpois(50, 2) | se(nope) ~ x, family = poisson(), data = d)))
say("se(s, sigma = TRUE) present", try_msg(class(frm(y | se(s, sigma = TRUE) ~ x, data = d))))
sd <- 0
say("se column named like a function (sd) missing from data, gaussian", try_msg(frm(y | se(sd) ~ x, data = d)))
say("student se missing col", try_msg(frm(y | se(nope) ~ x, family = student(), data = d)))
say("default_prior se column missing, weibull", try_msg(default_prior(abs(y) | se(nope) ~ x, family = weibull(), data = d)))

## 8. mixture(order =), family print, family(resp =) --------------------
hdr("8. mixture(order), print, family(resp)")
say("order none == NULL", identical(mixture(gaussian(), gaussian(), order = "none"), mixture(gaussian(), gaussian())))
say("order FALSE == NULL", identical(mixture(gaussian(), gaussian(), order = FALSE), mixture(gaussian(), gaussian())))
say("order mu", try_msg(mixture(gaussian(), gaussian(), order = "mu")))
say("order TRUE", try_msg(mixture(gaussian(), gaussian(), order = TRUE)))
say("order bad", try_msg(mixture(gaussian(), gaussian(), order = "x")))
say("order NA", try_msg(mixture(gaussian(), gaussian(), order = NA)))
say("order c()", try_msg(mixture(gaussian(), gaussian(), order = c("none", "mu"))))
say("nmix still refused", try_msg(mixture(gaussian(), nmix = 2)))
print(student(), links = TRUE)
print(mixture(gaussian(), exponential()), links = TRUE)
print(cumulative())
say("links bad", try_msg(print(student(), links = 3)))
set.seed(81)
d <- data.frame(y = rnorm(40), y2 = rpois(40, 2), x = rnorm(40))
fmv <- frm(bf(y ~ x) + bf(y2 ~ x, family = poisson()) + set_rescor(FALSE), data = d)
say("family(resp = 'y2')", try_msg(family(fmv, resp = "y2")$family))
say("family(resp = both) length", try_msg(length(family(fmv, resp = c("y", "y2")))))
say("family(resp = 'nope')", try_msg(family(fmv, resp = "nope")))
f1 <- frm(y ~ x, data = d)
say("univariate family(resp = 'y')", try_msg(family(f1, resp = "y")$family))

## 9. ranef(pars, groups) ------------------------------------------------
hdr("9. ranef(pars =, groups =)")
set.seed(91)
n <- 200
d <- data.frame(x = rnorm(n), g = factor(rep(1:10, 20)), h = factor(rep(1:5, 40)))
d$y <- rnorm(10)[d$g] + rnorm(5)[d$h] + d$x + rnorm(n)
f <- frm(bf(y ~ x + (1 + x | g) + (1 | h), sigma ~ (1 | g)), data = d)
r <- ranef(f)
say("names", names(r)); say("colnames g", colnames(r$g))
r2 <- ranef(f, pars = "x")
say("pars = 'x' names", names(r2)); say("pars = 'x' cols", lapply(r2, colnames))
r3 <- ranef(f, pars = "Intercept")
say("pars = 'Intercept' cols", lapply(r3, colnames))
r4 <- ranef(f, groups = "h")
say("groups = 'h'", names(r4))
say("groups = 'zz'", length(ranef(f, groups = "zz")))
r5 <- ranef(f, pars = "x", condVar = TRUE)
say("condSD dims follow", dim(attr(r5$g, "condSD")))
say("pars numeric", try_msg(ranef(f, pars = 1)))
say("values equal subset", isTRUE(all.equal(r2$g[, "x"], r$g[, "x"])))

## 10. data_name ----------------------------------------------------------
hdr("10. data_name")
set.seed(101)
dd <- data.frame(y = rnorm(40), x = rnorm(40))
f <- frm(y ~ x, data = dd)
say("name", attr(f$data, "data_name"))
f2 <- frm(y ~ x, data = dd[1:30, ])
say("call name", attr(f2$data, "data_name"))
new_data <- dd[5:40, ]
u <- update(f, newdata = new_data)
say("update(newdata = new_data)", attr(u$data, "data_name"))
u2 <- update(f, data = new_data)
say("update(data = new_data)", attr(u2$data, "data_name"))
f3 <- do.call(frm, list(y ~ x, data = dd))
say("do.call by value", attr(f3$data, "data_name"))
say("summary Data line", grep("Data", capture.output(summary(f3)), value = TRUE))
say("summary Data line named", grep("Data", capture.output(summary(f2)), value = TRUE))
g <- function(dat) frm(y ~ x, data = dat)
say("inside a function", attr(g(dd)$data, "data_name"))
long_name_for_a_data_frame_that_is_rather_long_indeed_yes <- dd
f4 <- frm(y ~ x, data = long_name_for_a_data_frame_that_is_rather_long_indeed_yes)
say("long name cut at 50", nchar(attr(f4$data, "data_name")))

## 11. default_prior on invalid responses; frm still validates -----------
hdr("11. default_prior / frm response validation")
set.seed(111)
d <- data.frame(x = rnorm(30), yb = rnorm(30), cnt = rnorm(30),
                y3 = sample(0:2, 30, TRUE), yo = sample(1:3, 30, TRUE))
for (fm in list(list(yb ~ x, Beta()), list(y3 ~ x, bernoulli()),
                list(cnt ~ x, poisson()), list(cnt ~ x, negbinomial()),
                list(cnt ~ x, lognormal()), list(yb ~ x, cumulative()),
                list(y3 | trials(1) ~ x, binomial()))) {
  cat("  ", deparse1(fm[[1]]), fm[[2]]$family, "\n")
  say("    default_prior rows", try_msg(nrow(default_prior(fm[[1]], data = d, family = fm[[2]]))))
  say("    frm", try_msg(class(frm(fm[[1]], data = d, family = fm[[2]]))))
  say("    validate_prior", try_msg(nrow(validate_prior(set_prior("normal(0,1)", class = "b"), fm[[1]], data = d, family = fm[[2]]))))
}

## 12. set_prior spellings reach the parameter ----------------------------
hdr("12. set_prior old and new spellings")
set.seed(121)
n <- 100
d <- data.frame(x = runif(n), z = rnorm(n))
d$y <- sin(3 * d$x) + 0.5 * d$z + rnorm(n, 0, 0.3)
base <- frm(y ~ z + s(x), data = d)
p_new <- frm(y ~ z + s(x), data = d, prior = set_prior("normal(5, 0.01)", class = "b", coef = "sx_1"))
p_old <- frm(y ~ z + s(x), data = d, prior = set_prior("normal(5, 0.01)", class = "b", coef = "s(x).fx1"))
say("fixef base / new / old", cbind(fixef(base)[, 1], fixef(p_new)[, 1], fixef(p_old)[, 1]))
say("old == new", identical(p_old$opt$par, p_new$opt$par))
tab <- default_prior(y ~ z + s(x), data = d)
tab$prior[tab$coef == "sx_1"] <- "normal(5, 0.01)"
p_tab <- try_msg(frm(y ~ z + s(x), data = d, prior = tab))
say("table edited on sx_1 row, round trip == new", if (is.character(p_tab)) p_tab else identical(p_tab$opt$par, p_new$opt$par))
vp <- try_msg(validate_prior(set_prior("normal(5, 0.01)", class = "b", coef = "s(x).fx1"), y ~ z + s(x), data = d))
say("validate_prior old spelling rows with prior", if (is.character(vp)) vp else as.data.frame(vp)[nzchar(vp$prior) & vp$prior != "(flat)", c("prior", "class", "coef")])
dn <- data.frame(x = rnorm(n)); dn$y <- 2 * exp(0.3 * dn$x) + rnorm(n, 0, 0.2)
nl_new <- try_msg(frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = dn,
                      prior = c(set_prior("normal(10, 0.01)", nlpar = "a", coef = "Intercept"),
                                set_prior("normal(0, 1)", nlpar = "b"))))
nl_old <- try_msg(frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = dn,
                      prior = c(set_prior("normal(10, 0.01)", nlpar = "a", coef = "(Intercept)"),
                                set_prior("normal(0, 1)", nlpar = "b"))))
say("nl Intercept new vs old", if (is.character(nl_new) || is.character(nl_old)) c(nl_new, nl_old) else c(fixef(nl_new)[1, 1], fixef(nl_old)[1, 1], identical(nl_new$opt$par, nl_old$opt$par)))
tabn <- default_prior(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = dn)
say("nl table", as.data.frame(tabn)[, c("prior", "class", "coef", "nlpar")])
# a numeric column named like x.fx1 (not a smooth)
d$v.fx1 <- rnorm(n)
tv <- try_msg(as.data.frame(default_prior(y ~ z + v.fx1, data = d))[, c("class", "coef")])
say("column named v.fx1 in prior table", tv)
pv <- try_msg(fixef(frm(y ~ z + v.fx1, data = d, prior = set_prior("normal(3, 0.01)", class = "b", coef = "v.fx1")))[, 1])
say("prior on v.fx1 applied", pv)

## 13. categorical numeric response --------------------------------------
hdr("13. categorical numeric response")
set.seed(131)
n <- 150
d <- data.frame(x = rnorm(n))
d$yn <- sample(c(2, 10, 5), n, TRUE)
d$yf <- factor(d$yn)
fn <- frm(yn ~ x, family = categorical(), data = d)
ff <- frm(yf ~ x, family = categorical(), data = d)
say("levels order", levels(ff$data$yf))
say("same objective", identical(fn$opt$objective, ff$opt$objective))
say("fitted dimnames", dimnames(fitted(fn))[[3]])
say("dpars", names(fn$estimates))
say("predict works", dim(try_msg(predict(fn, ndraws = 20))))
say("newdata fitted", dim(try_msg(fitted(fn, newdata = d[1:3, ]))))

## 14. gp summary: several gp terms, by-gp, approximate ------------------
hdr("14. gp summary")
set.seed(141)
n <- 80
d <- data.frame(x = runif(n, 0, 10), z = runif(n, 0, 5), f = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sin(d$x) + cos(d$z) + rnorm(n, 0, 0.3)
for (fo in list(y ~ gp(x), y ~ gp(x) + gp(z), y ~ gp(x, k = 10), y ~ gp(x, z), y ~ gp(x, by = f),
                bf(y ~ 1, sigma ~ gp(x)))) {
  cat("  ", deparse1(fo), "\n")
  ft <- try_msg(frm(fo, data = d))
  if (is.character(ft)) { say("   fit", ft); next }
  s <- try_msg(summary(ft)$gp)
  print(s)
}
