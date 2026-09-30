# Probe: can brms 2.23.0 sample with algorithm = "fixed_param" at a
# given parameter vector, and what does its smooth data look like.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
set.seed(1)
n <- 150
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sin(2 * pi * d$x) + d$z * (d$f == "b") + rnorm(n, 0, 0.3)
form <- brms::bf(y ~ f + s(x) + s(z, by = f))
sd <- standata(form, data = d)
print(names(sd))
print(colnames(sd$Xs)); print(dim(sd$Zs_1_1)); print(dim(sd$Zs_2_1))
print(colnames(sd$X))
cat(grep("parameters", stancode(form, data = d), value = TRUE), "\n")
sc <- stancode(form, data = d)
cat(sc)
