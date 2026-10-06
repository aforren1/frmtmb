# What brms 2.23.0 does when newdata holds a by-level of gp(x, by = f)
# that the fit never saw, measured on its own stored example fit6
# (volume ~ Trt + gp(Age, by = Trt)), which needs no compile.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
fit6 <- brms:::rename_pars(brms:::brmsfit_example6)
print(fit6$formula)
str(fit6$data)
print(head(variables(fit6), 40))
nd <- fit6$data[1:4, ]
nd$Trt <- factor(c("0", "1", "0", "1"), levels = levels(fit6$data$Trt))
ep <- posterior_epred(fit6, newdata = nd, resp = "volume")
cat("seen levels: dim", dim(ep), "\n")
nd2 <- nd
nd2$Trt <- factor(c("0", "1", "2", "1"))
r <- try(posterior_epred(fit6, newdata = nd2, resp = "volume"))
cat("new level as new factor level:\n"); print(r)
r2 <- try(posterior_epred(fit6, newdata = nd2, resp = "volume",
                          allow_new_levels = TRUE))
cat("with allow_new_levels:\n"); print(r2)
nd3 <- nd
nd3$Trt <- NULL
r3 <- try(posterior_epred(fit6, newdata = nd3, resp = "volume"))
cat("by column absent:\n"); print(class(r3))
# the posterior predictive at NEW Age positions: is the conditional
# GP draw random given the hyperparameters? Same draw twice, two seeds.
nd4 <- data.frame(Age = c(-1.9, 0.33, 3.1), Trt = factor(c("0", "0", "1"),
  levels = levels(fit6$data$Trt)), count = 0, visit = 1, patient = 1)
set.seed(1); a <- posterior_epred(fit6, newdata = nd4, resp = "volume",
                                  draw_ids = 1:5, re_formula = NA)
set.seed(2); b <- posterior_epred(fit6, newdata = nd4, resp = "volume",
                                  draw_ids = 1:5, re_formula = NA)
print(a); print(b)
cat("max |a - b|:", max(abs(a - b)), "\n")
