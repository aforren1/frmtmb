# Constructions behind the fit-requiring verdicts of the brms suite port,
# on the fixtures of tests/testthat/helper-brms-suite.R.
#   Rscript dev/brmsport-probe4.R > dev/brmsport-log/probe4.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/brmsport-lib", "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb); library(frmtmb.sample)})
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
show <- function(label, expr) {
  cat("\n--", label, "\n")
  r <- tryCatch(withCallingHandlers(expr,
         warning = function(w) {cat("  [warning]", conditionMessage(w), "\n"); invokeRestart("muffleWarning")},
         message = function(m) invokeRestart("muffleMessage")),
       error = function(e) {cat("  ERROR:", conditionMessage(e), "\n"); NULL})
  if (!is.null(r)) print(utils::head(r, 12))
  invisible(r)
}
fit1 <- brms_fixture(1); fit2 <- brms_fixture(2); fit3 <- brms_fixture(3)
fit4 <- brms_fixture(4); fit5 <- brms_fixture(5); fit6 <- brms_fixture(6)
show("(a) names(fit1)", names(fit1))
show("(a) fit1$data", fit1$data)
show("(a) class(model.frame(fit1)), dim", dim(model.frame(fit1)))
show("(b) hypothesis(fit3, 'Age') no relation", hypothesis(fit3, "Age")$hypothesis)
show("(b) hypothesis(fit3, 'b_Age x 0')", hypothesis(fit3, "b_Age x 0"))
show("(c) vcov(fit1) rownames", rownames(vcov(fit1)))
show("(d) names(summary(fit1))", names(summary(fit1)))
show("(d) class(summary(fit1))", class(summary(fit1)))
cat(utils::capture.output(print(summary(fit6)))[1:25], sep = "\n")
show("(e) residuals(fit4)", residuals(fit4))
show("(e) class pp_check(fit4, 'error_binned')", class(pp_check(fit4, "error_binned")))
show("(f) conditional_effects(fit1, c('Trtc', 'Trt'))", class(conditional_effects(fit1, effects = c("Trtc", "Trt"))))
ds1 <- brms_fixture_draws(1)
show("(g) formals posterior_epred.frmtmb_draws", names(formals(getS3method("posterior_epred", "frmtmb_draws"))))
show("(g) posterior_epred(ds1, point_estimate = 'median', ndraws_point_estimate = 2) dim",
     dim(posterior_epred(ds1, point_estimate = "median", ndraws_point_estimate = 2)))
show("(g) posterior_epred(ds1, not_an_argument = 2) dim", dim(posterior_epred(ds1, not_an_argument = 2)))
new_data <- data.frame(Age = rnorm(18), visit = rep(c(3, 2, 4), 6), Trt = rep(0:1, 9),
                       count = rep(c(5, 17, 28), 6), patient = rep(1:6, each = 3), Exp = 4, volume = 0)
show("(h) update(fit3, data = new_data3) nobs", {
  nd3 <- brms_fixture_data(3)[1:20, ]
  up <- update(fit3, data = nd3); c(nobs(fit3), nobs(up))})
show("(i) emmeans fit2 nlpar a", summary(emmeans::emmeans(fit2, "Age", nlpar = "a")))
show("(i) emmeans fit2 plain", summary(emmeans::emmeans(fit2, "Age")))
show("(i) emmeans fit6 epred", summary(emmeans::emmeans(fit6, "Age", by = "Trt", epred = TRUE)))
nd <- model.frame(fit1)[1:5, ]
show("(j) fitted(fit1, newdata = mf[1:5,])", fitted(fit1, newdata = nd))
nd$visit <- factor(c(1:4, 5))
show("(j) new visit level 5", fitted(fit1, newdata = nd))
show("(k) update(fit1, ~ . - Age + factor(Age))", class(update(fit1, ~ . - Age + factor(Age))))
show("(k) traceback-like: update(fit3, bf(~., family = acat()))", class(update(fit3, bf(~ ., family = acat()))))
nd2 <- data.frame(Age = c(0, -0.2), visit = c(1, 4), Trt = c(0, 1), count = c(20, 13),
                  patient = c(1, 42), Exp = c(2, 4), volume = 0)
show("(l) fitted newdata numeric Trt", fitted(fit1, newdata = nd2))
nd2$Trt <- factor(nd2$Trt); nd2$Exp <- factor(nd2$Exp, levels = 1:5, ordered = TRUE)
show("(l) fitted newdata factor Trt", fitted(fit1, newdata = nd2))
show("(m) str(fixef(fit1))", str(fixef(fit1)))
show("(m) str(ranef(fit1))", str(unclass(ranef(fit1)), max.level = 2))
show("(n) ngrps(fit1)", ngrps(fit1))
cat(utils::capture.output(print(fit1))[1:30], sep = "\n")
show("(p) pp_check ribbon_grouped", class(pp_check(fit1, "ribbon_grouped", group = "visit", x = "Age")))
show("(q) variables(fit4)", variables(fit4))
show("(q) coef names fit4", names(coef(fit4)))
show("(r) print(family(fit5))", print(family(fit5)))
show("(s) fitted(fit1) head", fitted(fit1))
