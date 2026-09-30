# A multivariate model with two ordinal responses, each with a threshold
# structure, and disc on one: names and the log-likelihood as the sum of
# the univariate fits' (no shared effect, so an identity). Also brms's
# own Stan code for the same model, which declares `delta` twice.
# Seed 20260930. Output: dev/ordinal-log-mv.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20260930)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
v <- stats::rlogis(n) - 0.5 * d$x
d$y2 <- 1L + (v > -1) + (v > 0.5)
f1 <- bf(y ~ x, family = cumulative(threshold = "equidistant"))
f2 <- bf(y2 ~ x, disc ~ 0 + z, family = sratio(threshold = "equidistant"))
fit <- frm(f1 + f2 + set_rescor(FALSE), data = d)
print(variables(fit))
u1 <- frm(y ~ x, family = cumulative(threshold = "equidistant"), data = d)
u2 <- frm(bf(y2 ~ x, disc ~ 0 + z),
          family = sratio(threshold = "equidistant"), data = d)
ll <- as.numeric(logLik(fit))
ls <- as.numeric(logLik(u1)) + as.numeric(logLik(u2))
cat("mv logLik", format(ll, digits = 12), " sum of univariate",
    format(ls, digits = 12), " relative difference",
    format(abs(ll - ls) / abs(ls), digits = 3), "\n")
print(summary(fit)$spec_pars)
cat("\n== brms 2.23.0 stancode, parameters block ==\n")
sc <- brms::stancode(brms::bf(y ~ x,
                              family = brms::cumulative(threshold =
                                                          "equidistant")) +
                       brms::bf(y2 ~ x, disc ~ 0 + z,
                                family = brms::sratio(threshold =
                                                        "equidistant")) +
                       brms::set_rescor(FALSE), data = d)
ln <- strsplit(sc, "\n")[[1]]
cat(grep("delta", ln, value = TRUE), sep = "\n")
r <- tryCatch({
  rstan::stanc(model_code = sc)
  "compiles"
}, error = function(e) paste("stanc ERROR:", conditionMessage(e)))
cat(substr(r, 1, 400), "\n")
