# Reviewer of lane defects: sections 4 to 14 of dev/defects-rev-probe.R,
# which stopped at section 3 on a probe error. Behavioral probes on the lane
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
say("names", names(r)); say("colnames", lapply(unclass(r), colnames))
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
say("values equal subset", isTRUE(all.equal(unclass(r2)[[1]][, "x"], unclass(r)[[1]][, "x"])))

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
