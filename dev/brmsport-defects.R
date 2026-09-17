# Constructions for the defects the brms suite port found, frmtmb first,
# brms 2.23.0 second on the same call. Sections S1 to S8 are the silent
# defects of dev/brmsport-findings.md section 5 in the same numbering;
# L1 to L8 are loud ones. brms runs on its own stored example fits
# (brms:::brmsfit_example1..6) or through make_standata(), so nothing
# compiles.
#
#   Rscript dev/brmsport-defects.R > dev/brmsport-log/defects.txt 2>&1
#
# Data seed 20260917 wherever data are simulated; the draws fixtures use
# frm_sample(seed = 20260917).
.libPaths(c("C:/Users/adf44/source/r/brmsport-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.sample)
})
h <- new.env()
sys.source("tests/testthat/helper-brms-suite.R", envir = h)
sys.source(file.path("extensions/frmtmb.sample/tests/testthat",
                     "helper-brms-suite-draws.R"), envir = h)
cat("frmtmb", format(packageVersion("frmtmb")), "frmtmb.sample",
    format(packageVersion("frmtmb.sample")), "brms",
    format(packageVersion("brms")), "\n")

show <- function(label, expr) {
  warns <- character()
  out <- tryCatch({
    v <- withCallingHandlers(expr, warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
    paste(utils::capture.output(print(v)), collapse = " | ")
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("  %-46s %s\n", label, substr(gsub(" +", " ", out), 1, 200)))
  if (length(warns)) {
    cat(sprintf("  %-46s [warning] %s\n", "",
                substr(unique(warns), 1, 150)))
  }
}
bex <- function(k) {
  brms:::rename_pars(get(paste0("brmsfit_example", k),
                         envir = asNamespace("brms")))
}
fit1 <- h$brms_fixture(1)
fit3 <- h$brms_fixture(3)
fit4 <- h$brms_fixture(4)
ds1 <- h$brms_fixture_draws(1)

cat("\n== S1 point_estimate and ndraws_point_estimate ignored on draws\n")
show("frmtmb.sample point_estimate = median, 2",
     dim(posterior_epred(ds1, point_estimate = "median",
                         ndraws_point_estimate = 2)))
show("brms point_estimate = median, 2",
     dim(brms::posterior_epred(bex(1), point_estimate = "median",
                               ndraws_point_estimate = 2)))
cat("  control: an unknown argument is ignored by BOTH (brms parity)\n")
show("frmtmb.sample not_an_argument = 2",
     dim(posterior_epred(ds1, not_an_argument = 2)))
show("brms not_an_argument = 2",
     dim(brms::posterior_epred(bex(1), not_an_argument = 2)))

cat("\n== S2 ar() and ma() with an expression as the time term\n")
set.seed(20260917)
d <- data.frame(g1 = rep(c(1, 1, 2, 2), each = 6),
                g2 = rep(c(1, 2, 1, 2), each = 6),
                t = c(1:6, 1:6, 1:6, 7:12))
d$x <- rnorm(24)
d$g <- interaction(d$g1, d$g2)
d$y <- d$x + as.numeric(stats::filter(rnorm(24), 0.6, "recursive"))
show("frmtmb ar(t, g, cov = TRUE) control",
     round(logLik(frm(y ~ x + ar(t, g, cov = TRUE), d)), 6))
show("frmtmb ar(x + t, g, cov = TRUE)",
     round(logLik(frm(y ~ x + ar(x + t, g, cov = TRUE), d)), 6))
show("frmtmb ar(t - 10 * x, g, cov = TRUE)",
     round(logLik(frm(y ~ x + ar(t - 10 * x, g, cov = TRUE), d)), 6))
show("frmtmb ma(t, g, cov = TRUE) control",
     round(logLik(frm(y ~ x + ma(t, g, cov = TRUE), d)), 6))
show("frmtmb ma(x + t, g, cov = TRUE)",
     round(logLik(frm(y ~ x + ma(x + t, g, cov = TRUE), d)), 6))
show("brms ar(x + t, g)",
     names(brms::make_standata(y ~ x + ar(x + t, g, cov = TRUE), d)))
show("brms ar(t - 10 * x, g)",
     names(brms::make_standata(y ~ x + ar(t - 10 * x, g, cov = TRUE), d)))
show("brms ma(x + t, g)",
     names(brms::make_standata(y ~ x + ma(x + t, g, cov = TRUE), d)))

cat("\n== S3 ar(gr = g1/g2) on numeric codes groups by the quotient\n")
# Series (2, 2) runs at times 7..12, so g1/g2 read as a division puts it
# in series (1, 1)'s group (both ratios are 1) with no repeated time.
show("frmtmb ar(t, gr = g1/g2, cov = TRUE)",
     round(logLik(frm(y ~ x + ar(t, gr = g1/g2, cov = TRUE), d)), 6))
show("frmtmb ar(t, gr = g1:g2, cov = TRUE)",
     round(logLik(frm(y ~ x + ar(t, gr = g1:g2, cov = TRUE), d)), 6))
show("brms ar(t, gr = g1/g2)",
     names(brms::make_standata(y ~ x + ar(t, gr = g1/g2, cov = TRUE), d)))

cat("\n== S4 fitted() and residuals() on draws\n")
show("frmtmb.sample fitted(ds1)", fitted(ds1))
show("frmtmb.sample residuals(ds1)", residuals(ds1))
show("brms dim(fitted(ex1))", dim(fitted(bex(1))))
show("brms dim(residuals(ex1))", dim(residuals(bex(1))))

cat("\n== S5 hypothesis() with no relation\n")
show("frmtmb 'Trt1 + Age' (silent)",
     hypothesis(fit3, "Trt1 + Age")$hypothesis)
show("frmtmb 'Trt1 + Age = 0'",
     hypothesis(fit3, "Trt1 + Age = 0")$hypothesis)
show("frmtmb 'b_Age x 0' (loud)", hypothesis(fit3, "b_Age x 0"))
show("brms hypothesis(ex1, 'Age')", brms::hypothesis(bex(1), "Age"))

cat("\n== S6 variables() of an ordinal fit omits thresholds\n")
show("frmtmb variables(fit4), sratio + cs()", variables(fit4))
show("frmtmb names(fit4$opt$par)", names(fit4$opt$par))
set.seed(20260917)
dc <- data.frame(y = factor(sample(1:3, 60, TRUE), ordered = TRUE),
                 x = rnorm(60))
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  fo <- frm(y ~ x, dc, family = get(fam)())
  show(sprintf("frmtmb %s: variables, opt$par", fam),
       c(variables = length(variables(fo)), par = length(fo$opt$par)))
}
dk <- data.frame(y = factor(sample(c("a", "b", "c"), 60, TRUE)),
                 x = rnorm(60))
fk <- frm(y ~ x, dk, family = categorical())
show("frmtmb categorical control: variables, par",
     c(variables = length(variables(fk)), par = length(fk$opt$par)))
show("brms variables(ex4) b_ and bcs_",
     grep("^(b_|bcs_)", brms::variables(bex(4)), value = TRUE))

cat("\n== S7 predictive errors on an ordinal fit\n")
r4 <- residuals(fit4)
p4 <- fitted(fit4)
y4 <- as.numeric(model.frame(fit4)$rating)
score <- as.numeric(p4 %*% seq_len(ncol(p4)))
show("frmtmb residuals(fit4)[1:4]", r4[1:4])
show("y code - sum(k P(Y = k)), [1:4]", (y4 - score)[1:4])
show("max |residual - (y - score)| / max |residual|",
     max(abs(r4 - (y4 - score))) / max(abs(r4)))
show("frmtmb pp_check(fit4, 'error_binned')",
     class(pp_check(fit4, "error_binned"))[1])
show("brms residuals(ex4)", dim(residuals(bex(4))))
show("brms pp_check(ex4, 'error_binned')",
     class(brms::pp_check(bex(4), "error_binned")))

cat("\n== S8 fit$data is the $ partial match of fit$data2\n")
show("names(fit1)", names(fit1))
show("fit1[['data', exact = TRUE]]", fit1[["data", exact = TRUE]])
show("identical(fit1$data, fit1$data2)", identical(fit1$data, fit1$data2))
set.seed(20260917)
A <- diag(6)
A[A == 0] <- 0.3
dimnames(A) <- list(1:6, 1:6)
dd <- data.frame(g = factor(rep(1:6, each = 5)), x = rnorm(30))
dd$y <- dd$x + rnorm(6)[dd$g] + rnorm(30)
fa <- frm(y ~ x + (1 | gr(g, cov = A)), dd, data2 = list(A = A))
show("fit with data2: names(fit$data)", names(fa$data))
show("fit with data2: dim(fit$data$A)", dim(fa$data$A))
show("brms dim(ex1$data)", dim(bex(1)$data))

cat("\n== L1 categorical() on integer codes\n")
show("frmtmb categorical, y = rep(1:3, 4)",
     frm(y ~ 1, data.frame(y = rep(1:3, 4)), family = categorical(),
         dry_run = "frame")$n_obs)
show("frmtmb categorical, y = factor(rep(1:3, 4))",
     frm(y ~ 1, data.frame(y = factor(rep(1:3, 4))),
         family = categorical(), dry_run = "frame")$n_obs)
show("brms categorical, y = rep(1:3, 4)",
     brms::make_standata(y ~ 1, data.frame(y = rep(1:3, 4)),
                         family = brms::categorical())$ncat)

cat("\n== L2 y ~ .\n")
set.seed(20260917)
d3 <- data.frame(y = rnorm(10), x1 = rnorm(10), x2 = rnorm(10))
show("frmtmb frm(y ~ ., d3)", coef(frm(y ~ ., d3)))
show("brms make_standata(y ~ ., d3)",
     colnames(brms::make_standata(y ~ ., d3)$X))

cat("\n== L3 a factor covariate in a nonlinear body\n")
d2 <- h$brms_fixture_data(2)
show("frmtmb Trt as the factor",
     logLik(frm(bf(count ~ 1 / (1 + exp(-a)) * exp(b * Trt), a ~ Age,
                   b ~ Age, nl = TRUE), d2, family = Gamma("identity"))))
show("brms C_1 for the factor",
     head(brms::make_standata(
       brms::bf(count ~ 1 / (1 + exp(-a)) * exp(b * Trt), a ~ Age,
                b ~ Age, nl = TRUE), d2,
       family = brms::brmsfamily("gamma", "identity"))$C_1))

cat("\n== L4 newdata with a factor's numeric codes\n")
nd <- data.frame(Age = c(0, -0.2), visit = c(1, 4), Trt = c(0, 1),
                 count = c(20, 13), patient = c(1, 42), Exp = c(2, 4),
                 volume = 0)
show("frmtmb fitted(fit1, newdata = nd)", fitted(fit1, newdata = nd))
show("brms dim(fitted(ex1, newdata = nd))",
     dim(fitted(bex(1), newdata = nd)))

cat("\n== L5 fit-surface shapes (user rule 3, item 2.6f)\n")
show("frmtmb dim(fitted(fit1))", dim(fitted(fit1)))
show("brms dim(fitted(ex1))", dim(fitted(bex(1))))
show("frmtmb class(fixef(fit1))", class(fixef(fit1)))
show("brms dim(fixef(ex1))", dim(fixef(bex(1))))
show("frmtmb ngrps(fit1)", ngrps(fit1))
show("brms ngrps(ex1)", brms::ngrps(bex(1)))
show("frmtmb dim(vcov(fit1))", dim(vcov(fit1)))
show("brms dim(vcov(ex1))", dim(vcov(bex(1))))
show("frmtmb names(summary(fit1))", names(summary(fit1)))

cat("\n== L6 update(data =)\n")
show("frmtmb update(fit3, data = rows 1..20) nobs",
     nobs(update(fit3, data = h$brms_fixture_data(3)[1:20, ])))
show("brms update(ex1, data = ...)",
     update(bex(1), data = bex(1)$data[1:20, ], testmode = TRUE))

cat("\n== L7 emmeans on a mo() term\n")
show("frmtmb emmeans(fit1, 'Exp')",
     nrow(summary(emmeans::emmeans(fit1, "Exp"))))
show("brms emmeans(ex1, 'Exp')",
     nrow(summary(emmeans::emmeans(bex(1), "Exp"))))

cat("\n== L8 default_prior() validates the response (pre-existing)\n")
set.seed(20260917)
dp <- data.frame(y = rnorm(10), z = rnorm(10), g = rep(1:2, 5))
show("frmtmb default_prior Beta on rnorm y",
     nrow(default_prior(bf(y ~ 1, phi ~ z + (1 | g), family = Beta()), dp)))
show("brms default_prior Beta on rnorm y",
     nrow(brms::default_prior(brms::bf(y ~ 1, phi ~ z + (1 | g),
                                       family = brms::Beta()), dp)))
