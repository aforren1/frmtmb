# Constructions behind the no-fit verdicts of the brms suite port: for
# each questionable assertion, what frmtmb does and what brms 2.23.0
# does on the same call, so a verdict of defect or divergence rests on
# both answers. brms runs through make_standata()/default_prior()/
# brmsterms(), which never compile Stan.
#
#   Rscript dev/brmsport-probe3.R > dev/brmsport-log/probe3.txt 2>&1
lib <- "C:/Users/adf44/source/r/brmsport-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "brms",
    format(packageVersion("brms")), "\n")

show <- function(label, expr) {
  out <- tryCatch({
    v <- withCallingHandlers(expr, warning = function(w) {
      cat("  [warning]", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }, message = function(m) {
      cat("  [message]", conditionMessage(m))
      invokeRestart("muffleMessage")
    })
    paste("OK:", paste(utils::capture.output(str(v, max.level = 1,
                                                 give.attr = FALSE))[1:2],
                       collapse = " | "))
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 220)))
  cat(sprintf("%-58s %s\n", label, out))
}

set.seed(20260917)
dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2),
                  sei = rexp(10), g1 = rep(1:2, 5), g2 = rep(1:5, 2))
pdat <- transform(dat, yp = rexp(10) + 0.1, yc = rpois(10, 3))

cat("\n-- brm:81 se() with weibull, sei present\n")
show("frm weibull se(sei)",
     frm(yp | se(sei) ~ x, pdat, family = weibull(), dry_run = "frame"))
show("brms weibull se(sei)",
     brms::make_standata(yp | se(sei) ~ x, pdat, family = brms::weibull()))

cat("\n-- brm:106/108 ar() with cov = TRUE\n")
show("frm ar(x + y, g, cov = TRUE)",
     frm(y ~ ar(x + y, g, cov = TRUE), dat, dry_run = "frame"))
show("frm ar(gr = g1/g2, cov = TRUE)",
     frm(y ~ ar(gr = g1/g2, cov = TRUE), dat, dry_run = "frame"))
show("frm ar(time = x, gr = g, cov = TRUE) control",
     frm(y ~ ar(x, g, cov = TRUE), dat, dry_run = "frame"))

cat("\n-- brm:102 set_rescor() on the right-hand side\n")
show("frm y ~ 1 + set_rescor(TRUE)", frm(y ~ 1 + set_rescor(TRUE), dat))

cat("\n-- brm:116 and standata:83 categorical on integer codes\n")
inh <- brms::inhaler
show("frm categorical rating ~ treat",
     frm(rating ~ treat, inh, family = categorical(), dry_run = "frame"))
show("frm categorical factor(rating) ~ treat",
     frm(factor(rating) ~ treat, transform(inh, rating = factor(rating)),
         family = categorical(), dry_run = "frame"))
show("brms categorical rating ~ treat",
     brms::make_standata(rating ~ treat, inh, family = brms::categorical()))
show("frm categorical y = rep(1:10, 5)",
     frm(y ~ 1, data.frame(y = rep(1:10, 5)), family = categorical(),
         dry_run = "frame"))

cat("\n-- standata:75 bernoulli on two negative codes\n")
show("frm bernoulli y in {-1,-2}",
     frm(y ~ 1, data.frame(y = rep(-c(1:2), 5)), family = bernoulli(),
         dry_run = "frame"))
show("brms bernoulli y in {-1,-2}",
     brms::make_standata(y ~ 1, data.frame(y = rep(-c(1:2), 5)),
                         family = brms::bernoulli())$Y)

cat("\n-- standata:927 dot expansion\n")
d3 <- data.frame(y = rnorm(10), x1 = rnorm(10), x2 = rnorm(10))
show("frm y ~ . (fit)", coef(frm(y ~ ., d3)))
show("frm y ~ . (frame)", frm(y ~ ., d3, dry_run = "frame"))
show("brms y ~ .", colnames(brms::make_standata(y ~ ., d3)$X))

cat("\n-- standata:718 poly() column names\n")
d4 <- data.frame(y = rnorm(10), x = rnorm(10))
f4 <- frm(y ~ 1 + poly(x, 3), d4)
show("frm variables(fit)", variables(f4))
show("frm fixef names", rownames(fixef(f4)$mu))

cat("\n-- standata:724 fixed nu: where the frame keeps it\n")
fr <- frm(bf(y ~ 1, nu = 3), list(y = 1:10), family = student(),
          dry_run = "frame")
show("spec fixed", fr$spec$responses$y$fixed)
str(fr$spec$responses$y, max.level = 1)

cat("\n-- brmsterms:44/46 fixed dpar range, through frm()\n")
dg <- data.frame(y = rexp(10) + 0.1)
show("frm Gamma shape = -2", frm(bf(y ~ 1, shape = -2), dg, family = Gamma()))
show("brms Gamma shape = -2",
     brms::brmsterms(brms::bf(y ~ 1, shape = -2, family = Gamma())))
show("frm asym_laplace quantile = 1.5",
     frm(bf(y ~ 1, quantile = 1.5), dg, family = asym_laplace()))
show("frm asym_laplace quantile = 0.5 control",
     frm(bf(y ~ 1, quantile = 0.5), dg, family = asym_laplace()))
show("frm gaussian sigma = 4", frm(bf(y ~ 1, sigma = 4), dg))
show("frm zero_inflated_beta zi = 0.5",
     frm(bf(y ~ 1, zi = 0.5), data.frame(y = c(0, runif(9))),
         family = zero_inflated_beta(), dry_run = "frame"))

cat("\n-- brmsterms:89 unused = ~ x\n")
du <- data.frame(y = rnorm(10), x = rnorm(10))
show("frm bf(y ~ 1, unused = ~ x)", frm(bf(y ~ 1, unused = ~x), du))

cat("\n-- priors:67 default_prior on a Beta dpar model with rnorm y\n")
dp <- data.frame(y = rnorm(10), x = rnorm(10), z = rnorm(10),
                 g = rep(1:2, 5))
show("frm default_prior Beta",
     default_prior(bf(y ~ 1, phi ~ z + (1 | g), family = Beta()), data = dp))
show("brms default_prior Beta",
     brms::default_prior(brms::bf(y ~ 1, phi ~ z + (1 | g),
                                  family = brms::Beta()), data = dp))
show("frm get_prior Beta",
     get_prior(bf(y ~ 1, phi ~ z + (1 | g), family = Beta()), data = dp))

cat("\n-- standata:179 se() of zero\n")
show("frm se(s) with a zero",
     frm(y | se(s) ~ 1, data.frame(y = rnorm(3), s = c(0, 1, 2)),
         dry_run = "frame"))
show("brms se(s) with a zero",
     brms::make_standata(y | se(s) ~ 1,
                         data.frame(y = rnorm(3), s = c(0, 1, 2)))$se)

cat("\n-- standata:970 reserved Intercept\n")
show("frm y ~ 0 + Intercept",
     frm(y ~ 0 + Intercept, data.frame(y = 1:10 + rnorm(10))))
