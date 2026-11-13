# Lane names, backlog 3164: what brms 2.23.0 does with a categorical
# response holding an unused factor level, with a named multinomial
# response matrix, and how it names posterior_epred()'s margins on a
# categorical fit. Output: dev/names-out/brms-3164.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
sink("dev/names-out/brms-3164.txt", split = TRUE)
set.seed(3164)
n <- 90
d <- data.frame(x = rnorm(n))
d$yc <- factor(sample(c("a", "b", "c"), n, TRUE),
               levels = c("a", "b", "c", "unused"))
cat("levels(d$yc):", levels(d$yc), "\n")
sd1 <- standata(bf(yc ~ x), data = d, family = categorical())
cat("categorical with an unused level: ncat =", sd1$ncat, "\n")
cat("stancode dpars:",
    grep("^  vector\\[N\\] mu", strsplit(stancode(bf(yc ~ x), data = d,
                                                 family = categorical()),
                                        "\n")[[1]], value = TRUE), "\n")
d$tr <- 10
d$y1 <- rbinom(n, 10, 0.3)
d$y2 <- rbinom(n, 10 - d$y1, 0.5)
d$y3 <- d$tr - d$y1 - d$y2
d$Y <- with(d, cbind(y1, y2, y3))
sc <- strsplit(stancode(bf(Y | trials(tr) ~ x), data = d,
                        family = multinomial()), "\n")[[1]]
cat("multinomial named columns, dpars:\n")
cat(grep("vector\\[N\\] mu", sc, value = TRUE), sep = "\n")
fit <- suppressWarnings(suppressMessages(brm(
  bf(yc ~ x), data = droplevels(d), family = categorical(), chains = 1,
  iter = 60, warmup = 30, seed = 1, refresh = 0, backend = "rstan")))
pe <- posterior_epred(fit)
cat("posterior_epred dims:", dim(pe), "\n")
str(dimnames(pe))
cat("posterior_predict dimnames:\n")
str(dimnames(posterior_predict(fit)))
fit2 <- tryCatch(suppressWarnings(suppressMessages(brm(
  bf(yc ~ x), data = d, family = categorical(), chains = 1,
  iter = 60, warmup = 30, seed = 1, refresh = 0, backend = "rstan"))),
  error = function(e) e)
if (inherits(fit2, "error")) {
  cat("fit with unused level: ERROR", conditionMessage(fit2), "\n")
} else {
  cat("fit with unused level, variables:",
      grep("^b_", variables(fit2), value = TRUE), "\n")
}
sink()
