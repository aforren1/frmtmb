# Item 2: default effects of conditional_effects() on an offset model.
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(21)
n <- 80
d <- data.frame(x = rnorm(n), time = runif(n, 1, 3), f = gl(2, 40))
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x) * d$time)
f <- frm(bf(yc ~ x + f + offset(log(time))), data = d, family = poisson())
print(names(tryCatch(conditional_effects(f), error = function(e) conditionMessage(e))))
print(names(model.frame(f)))
f2 <- frm(bf(yc ~ x * f + offset(log(time))), data = d, family = poisson())
print(names(tryCatch(conditional_effects(f2), error = function(e) conditionMessage(e))))
