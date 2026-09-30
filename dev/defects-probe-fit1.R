# Lane wt-defects: does brms's own fixture-1 formula, arma() with its
# default cov = FALSE, fit and answer predict(newdata = )?
source("dev/defects-pre.R")
d <- brms_fixture_data(1)
f <- bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) +
          (1 + Trt | visit) + arma(visit, patient), sigma ~ Trt)
fit <- show(frm(f, data = d, family = student()))
print(fit$opt$convergence); print(fit$opt$message)
nd <- data.frame(Age = c(0, -0.2), visit = c(1, 4), Trt = c(1, 0),
                 count = c(2, 10), patient = c(1, 42), Exp = c(1, 2),
                 volume = 0)
show(print(dim(predict(fit, newdata = nd))))
nd$visit <- c(1, 6)
show(print(dim(predict(fit, newdata = nd, allow_new_levels = TRUE))))
df <- d[1:10, ]; df$count[8:10] <- NA
show(print(predict(fit, newdata = df, ndraws = 1)[, "Estimate"]))
