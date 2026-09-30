# Lane wt-defects: every defect row this lane owns, on brms's own call.
#   DEFECTS_ARM=before Rscript dev/defects-probe-all.R   # rellib-r3
#   Rscript dev/defects-probe-all.R                      # lane library
source("dev/defects-pre.R")
set.seed(20260929)
fit1 <- brms_fixture(1); fit2 <- brms_fixture(2); fit3 <- brms_fixture(3)
fit4 <- brms_fixture(4); fit5 <- brms_fixture(5); fit6 <- brms_fixture(6)
dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))
row <- function(id) cat("\n==", id, "\n")

row("brm:81"); show(frm(y | se(sei) ~ x, dat, family = weibull()))
row("brm:102"); show(frm(y ~ 1 + set_rescor(TRUE), data = dat))
row("brm:116")
show(frm(rating ~ treat + (cs(period) | subject), data = brms::inhaler,
         family = categorical()))
row("methods:283"); show(print(family(fit6, resp = "count")))
row("methods:284"); show(print(family(fit1), links = TRUE))
row("methods:285"); show(print(family(fit5)))
nd <- data.frame(Age = c(0, -0.2), visit = c(1, 4), Trt = c(0, 1),
                 count = c(20, 13), patient = c(1, 42), Exp = c(2, 4),
                 volume = 0)
row("methods:301"); show(print(dim(fitted(fit1, newdata = nd))))
nd$visit <- c(1, 6)
row("methods:305")
show(print(dim(fitted(fit1, newdata = nd, allow_new_levels = TRUE))))
nd2 <- data.frame(Age = 0, visit = paste0("a", 1:100), Trt = 0, count = 20,
                  patient = 1, Exp = 2, volume = 0)
row("methods:314")
show(print(dim(fitted(fit1, newdata = nd2, allow_new_levels = TRUE,
                      sample_new_levels = "old_levels", ndraws = 10))))
row("methods:317")
show(print(dim(fitted(fit1, newdata = nd2, allow_new_levels = TRUE,
                      sample_new_levels = "gaussian", ndraws = 1))))
row("methods:352")
show(print(dim(fitted(fit4, newdata = fit4$data[1, ], scale = "linear"))))
show(print(dim(fitted(fit4, newdata = fit4$data[1, ]))))
row("methods:391")
hyp <- hypothesis(fit1, c("Age > Trt1", "Trt1:Age = -1"))
show(print(class(plot(hyp, plot = FALSE)[[1]])))
row("methods:396")
hyp <- hypothesis(fit1, "Intercept = 0", class = "sd", group = "visit")
show(print(class(plot(hyp, ignore_prior = TRUE, plot = FALSE)[[1]])))
nd <- data.frame(Age = c(0, -0.2), visit = c(1, 4), Trt = c(1, 0),
                 count = c(2, 10), patient = c(1, 42), Exp = c(1, 2),
                 volume = 0)
row("methods:737"); show(print(dim(predict(fit1, newdata = nd))))
nd$visit <- c(1, 6)
row("methods:741")
show(print(dim(predict(fit1, newdata = nd, allow_new_levels = TRUE))))
row("methods:747")
df <- fit1$data[1:10, ]; df$count[8:10] <- NA
show(print(predict(fit1, newdata = df, ndraws = 1)[, "Estimate"]))
row("methods:772")
nd5 <- fit5$data[1:5, ]; nd5$patient <- "a"
show(print(dim(predict(fit5, nd5, allow_new_levels = TRUE,
                       sample_new_levels = "old_levels"))))
row("methods:814"); show(print(dimnames(ranef(fit1, pars = "Trt1")$visit)))
row("methods:817"); show(print(length(ranef(fit1, groups = "a"))))
epilepsy <- brms::epilepsy
nd <- cbind(epilepsy[1:10, ], Exp = rep(1:5, 2), volume = 0)
row("methods:828"); show(print(dim(residuals(fit1, newdata = nd))))
nd$visit <- rep(1:5, 2)
row("methods:832")
show(print(dim(residuals(fit1, newdata = nd, allow_new_levels = TRUE))))
row("methods:840"); show(print(dim(residuals(fit6))))
row("methods:904"); show(print(summary(fit6)))
new_data <- data.frame(Age = rnorm(18), visit = rep(c(3, 2, 4), 6),
                       Trt = rep(0:1, 9), count = rep(c(5, 17, 28), 6),
                       patient = rep(1:6, each = 3), Exp = 4, volume = 0)
row("methods:924")
up <- show(update(fit1, newdata = new_data))
show(print(attr(up$data, "data_name")))
row("methods:927"); show(print(class(update(fit1, data = new_data))))
row("methods:995"); show(print(identical(variables(fit3), parnames(fit3))))
row("data-helpers:7")
fit <- fit1; newdata <- fit$data[1:5, ]
show(invisible(fitted(fit, newdata = newdata)))
row("families:102"); show(mixture(poisson, binomial, order = "x"))
d10 <- data.frame(y = rnorm(10), x = rnorm(10), z = rnorm(10),
                  g = rep(1:2, 5))
row("priors:55")
p <- show(default_prior(y ~ z + s(x) + (1 | g), data = d10))
show(print(p[p$class == "b", ]$coef))
row("priors:74")
show(print(default_prior(bf(y ~ 1, phi ~ z + (1 | g), family = Beta()),
                         data = d10)))
row("priors:101")
dc <- data.frame(y2 = c(1, rep(1:3, 3)), x = rnorm(10), g = rep(1:2, 5))
p <- show(default_prior(y2 ~ x + (x | ID1 | g), data = dc,
                        family = categorical()))
show(print(p[p$dpar == "mu2" & p$class == "b", "coef"]))
row("standata:75")
show(print(frm(y ~ 1, data = data.frame(y = rep(-c(1:2), 5)),
               family = bernoulli(), dry_run = "frame")$y))
row("standata:83")
show(print(frm(y ~ 1, data = data.frame(y = rep(11:20, 5)),
               family = categorical(), dry_run = "frame")$y))
row("standata:142")
show(frm(y ~ 1, data = data.frame(y = rep(0:1, 5)),
         family = categorical(), dry_run = "frame"))
