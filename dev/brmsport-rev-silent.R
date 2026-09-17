# Reviewer (brmsport): the eight silent claims on my own constructions,
# frmtmb first and brms 2.23.0 second, plus candidates the lane missed.
#   Rscript dev/brmsport-rev-silent.R > dev/brmsport-log/rev-silent.txt 2>&1
# Data seed 91 wherever data are simulated.
.libPaths(c("C:/Users/adf44/source/r/brmsport-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat); library(frmtmb); library(frmtmb.sample)
})
h <- new.env()
sys.source("tests/testthat/helper-brms-suite.R", envir = h)
sys.source(file.path("extensions/frmtmb.sample/tests/testthat",
                     "helper-brms-suite-draws.R"), envir = h)
show <- function(label, expr) {
  warns <- character()
  out <- tryCatch({
    v <- withCallingHandlers(expr, warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
    paste(utils::capture.output(print(v)), collapse = " | ")
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("  %-50s %s\n", label, substr(gsub(" +", " ", out), 1, 300)))
  if (length(warns)) cat(sprintf("  %-50s [warning] %s\n", "",
                                 substr(unique(warns), 1, 160)))
}
bex <- function(k) brms:::rename_pars(get(paste0("brmsfit_example", k),
                                          envir = asNamespace("brms")))

cat("\n== R1 dots on draws methods: does brms swallow too?\n")
ds1 <- h$brms_fixture_draws(1)
b1 <- bex(1)
for (m in c("posterior_epred", "posterior_linpred", "posterior_predict",
            "predictive_error", "log_lik")) {
  show(paste("frmtmb.sample", m, "(ndraw = 5) typo"),
       dim(get(m)(ds1, ndraw = 5)))
  show(paste("brms", m, "(ndraw = 5) typo"),
       dim(get(m, envir = asNamespace("brms"))(b1, ndraw = 5)))
  show(paste("brms", m, "(not_an_argument = 2)"),
       dim(get(m, envir = asNamespace("brms"))(b1, not_an_argument = 2)))
}
show("frmtmb.sample posterior_epred(ds, ndraws = 5)",
     dim(posterior_epred(ds1, ndraws = 5)))
show("brms nrow epred(point_estimate='median', 2)",
     dim(brms::posterior_epred(b1, point_estimate = "median",
                               ndraws_point_estimate = 2)))

cat("\n== R2/R3 ar() term arguments\n")
set.seed(91)
d <- data.frame(g1 = rep(c(1, 1, 2, 2), each = 8),
                g2 = rep(c(1, 2, 1, 2), each = 8),
                t = rep(1:8, 4))
d$t[25:32] <- 9:16
d$x <- rnorm(32)
d$g <- interaction(d$g1, d$g2)
d$y <- d$x + as.numeric(stats::filter(rnorm(32), 0.5, "recursive"))
show("frmtmb ar(t, g) control logLik",
     logLik(frm(y ~ x + ar(t, g, cov = TRUE), d)))
show("frmtmb ar(t, g2) control (merged groups) logLik",
     logLik(frm(y ~ x + ar(t, g2, cov = TRUE), d)))
show("frmtmb ar(2 * t, g) logLik",
     logLik(frm(y ~ x + ar(2 * t, g, cov = TRUE), d)))
show("frmtmb ar(t - 10 * x, g) logLik",
     logLik(frm(y ~ x + ar(t - 10 * x, g, cov = TRUE), d)))
show("frmtmb ma(x + t, g) logLik",
     logLik(frm(y ~ x + ma(x + t, g, cov = TRUE), d)))
show("frmtmb ar(t, gr = g1 * g2) logLik",
     logLik(frm(y ~ x + ar(t, gr = g1 * g2, cov = TRUE), d)))
show("frmtmb ar(t, gr = g1 + g2) logLik",
     logLik(frm(y ~ x + ar(t, gr = g1 + g2, cov = TRUE), d)))
show("frmtmb ar(t, gr = g1/g2) logLik",
     logLik(frm(y ~ x + ar(t, gr = g1/g2, cov = TRUE), d)))
show("brms ar(t - 10 * x, g)",
     names(brms::make_standata(y ~ x + ar(t - 10 * x, g, cov = TRUE), d)))
show("brms ar(t, gr = g1 + g2)",
     names(brms::make_standata(y ~ x + ar(t, gr = g1 + g2, cov = TRUE), d)))
show("brms ar(t, gr = g1 * g2)",
     names(brms::make_standata(y ~ x + ar(t, gr = g1 * g2, cov = TRUE), d)))
show("brms ma(x + t, g)",
     names(brms::make_standata(y ~ x + ma(x + t, g, cov = TRUE), d)))

cat("\n== R4 fitted()/residuals() on draws\n")
show("frmtmb.sample fitted(ds1)", fitted(ds1))
show("frmtmb.sample residuals(ds1)", residuals(ds1))
show("frmtmb.sample class(ds1)", class(ds1))
show("brms dim(residuals(ex1))", dim(residuals(b1)))

cat("\n== R5 hypothesis() with no relation\n")
fit1 <- h$brms_fixture(1); fit3 <- h$brms_fixture(3)
show("frmtmb hypothesis(fit3, 'Trt1 + Age')",
     hypothesis(fit3, "Trt1 + Age")$hypothesis)
show("frmtmb hypothesis(fit3, 'Trt1 + Age = 0')",
     hypothesis(fit3, "Trt1 + Age = 0")$hypothesis)
show("frmtmb hypothesis(fit3, 'Age x 0')",
     hypothesis(fit3, "Age x 0")$hypothesis)
show("frmtmb hypothesis(fit3, 'Age 0')", hypothesis(fit3, "Age 0")$hypothesis)
show("brms hypothesis(ex3, 'Trt1 + Age')",
     brms::hypothesis(bex(3), "Trt1 + Age"))

cat("\n== R6 variables() completeness\n")
fit4 <- h$brms_fixture(4); fit5 <- h$brms_fixture(5)
show("frmtmb variables(fit1)", variables(fit1))
want1 <- c("b_Intercept", "bsp_moExp", "ar[1]", "cor_visit__Intercept__Trt1",
           "nu")
show("brms fit-side names of 978 absent from variables(fit1)",
     setdiff(want1, variables(fit1)))
fit2 <- h$brms_fixture(2)
show("frmtmb variables(fit2)", variables(fit2))
want2 <- c("b_a_Intercept", "b_b_Age", "sd_patient__b_Intercept",
           "cor_patient__a_Intercept__b_Intercept")
show("fit-side names of 984 absent from variables(fit2)",
     setdiff(want2, variables(fit2)))
show("frmtmb variables(fit5)", variables(fit5))
show("frmtmb length(fit5$opt$par)", length(fit5$opt$par))
show("frmtmb variables(fit4)", variables(fit4))
set.seed(91)
dc <- data.frame(y = factor(sample(1:3, 60, TRUE), ordered = TRUE),
                 x = rnorm(60))
fc <- frm(y ~ x, dc, family = cumulative())
show("frmtmb variables(cumulative fit)", variables(fc))
show("frmtmb names(opt$par) cumulative", names(fc$opt$par))
dk <- data.frame(y = factor(sample(c("a", "b", "c"), 60, TRUE)),
                 x = rnorm(60))
fk <- frm(y ~ x, dk, family = categorical())
show("frmtmb variables(categorical fit)", variables(fk))
show("frmtmb names(opt$par) categorical", names(fk$opt$par))
show("brms variables(ex4)", brms::variables(bex(4)))

cat("\n== R7 residuals() on an ordinal fit: what is the number\n")
r4 <- residuals(fit4)
mf4 <- model.frame(fit4)
show("frmtmb head residuals(fit4)", head(r4))
show("frmtmb head y codes", head(as.integer(mf4[[1]])))
show("frmtmb head fitted(fit4)", head(fitted(fit4)))
show("y - sum(k * P(k))", head(as.integer(mf4[[1]]) -
                               as.vector(fitted(fit4) %*% 1:4)))

cat("\n== R8 fit$data mechanism\n")
show("names(fit1)", names(fit1))
show("exists $.frmtmb_fit", exists("$.frmtmb_fit", envir = asNamespace("frmtmb")))
show("fit1[['data', exact = TRUE]]", fit1[["data", exact = TRUE]])
set.seed(91)
ng <- 6
A <- diag(ng); A[A == 0] <- 0.3; dimnames(A) <- list(1:ng, 1:ng)
dd <- data.frame(g = factor(rep(1:ng, each = 5)), x = rnorm(30))
dd$y <- dd$x + rnorm(ng)[dd$g] + rnorm(30)
fa <- frm(y ~ x + (1 | gr(g, cov = A)), dd, data2 = list(A = A))
show("fit with data2: fit$data", str(fa$data))
show("fit with data2: nrow(fit$data)", nrow(fa$data))
show("fit with data2: names(fit$data)", names(fa$data))

cat("\n== R9 hollow-pass candidates\n")
hyp <- hypothesis(fit1, "Intercept = 0", class = "sd", group = "visit")
show("brmsfit-methods:394 Evid.Ratio[1]", hyp$hypothesis$Evid.Ratio[1])
show("  is.numeric", is.numeric(hyp$hypothesis$Evid.Ratio[1]))
show("brms Evid.Ratio[1] on ex1",
     brms::hypothesis(b1, "Intercept = 0", class = "sd",
                      group = "visit")$hypothesis$Evid.Ratio[1])
fit2 <- h$brms_fixture(2)
show("brmsfit-methods:698 pp_check(fit2, 'violin_grouped')",
     pp_check(fit2, "violin_grouped"))
show("brms pp_check(ex2, 'violin_grouped')",
     brms::pp_check(bex(2), "violin_grouped"))
show("fitted(fit1, newdata = mf[1:5,]) vs fitted(fit1)[1:5]", {
  a <- fitted(fit1, newdata = model.frame(fit1)[1:5, ])
  b <- fitted(fit1)[1:5]
  c(max_abs_diff = max(abs(a - b)), rel = max(abs(a - b) / abs(b)))
})
