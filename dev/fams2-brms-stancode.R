# Generated Stan code and data for the families this lane adds (brms 2.23.0).
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
n <- 40
d <- data.frame(x = rnorm(n), g = gl(4, 10))
d$yx <- pmin(1, pmax(0, runif(n, -0.1, 1.1)))
d$tr <- 10L
d$ybb <- rbinom(n, 10, 0.3)
d$yhc <- sample(0:4, n, TRUE)
show <- function(lab, expr) {
  cat("\n=====", lab, "=====\n")
  print(expr)
}
show("xbeta", stancode(bf(yx ~ x), data = d, family = xbeta()))
show("xbeta kappa ~ x", stancode(bf(yx ~ x, kappa ~ x), data = d,
                                 family = xbeta()))
show("zibb", stancode(bf(ybb | trials(tr) ~ x), data = d,
                      family = zero_inflated_beta_binomial()))
show("zibb zi~x", stancode(bf(ybb | trials(tr) ~ x, zi ~ x), data = d,
                      family = zero_inflated_beta_binomial()))
show("bb", stancode(bf(ybb | trials(tr) ~ x), data = d,
                    family = beta_binomial()))
show("hurdle_cumulative", stancode(bf(yhc ~ x), data = d,
                                   family = hurdle_cumulative()))
show("hurdle_cumulative hu ~ x", stancode(bf(yhc ~ x, hu ~ x), data = d,
                                   family = hurdle_cumulative()))
sd <- standata(bf(yhc ~ x), data = d, family = hurdle_cumulative())
show("hurdle_cumulative standata Y / nthres", list(Y = sd$Y, nthres = sd$nthres))
show("hurdle_cumulative default_prior",
     default_prior(bf(yhc ~ x), data = d, family = hurdle_cumulative()))
show("xbeta default_prior",
     default_prior(bf(yx ~ x), data = d, family = xbeta()))
show("zibb default_prior",
     default_prior(bf(ybb | trials(tr) ~ x), data = d,
                   family = zero_inflated_beta_binomial()))
d$yo <- sample(1:4, n, TRUE)
show("sratio equidistant", stancode(bf(yo ~ x), data = d,
                                    family = sratio(threshold = "equidistant")))
show("cumulative equidistant", stancode(bf(yo ~ x), data = d,
                                    family = cumulative(threshold = "equidistant")))
show("cumulative sum_to_zero", stancode(bf(yo ~ x), data = d,
                                    family = cumulative(threshold = "sum_to_zero")))
show("acat sum_to_zero", stancode(bf(yo ~ x), data = d,
                                    family = acat(threshold = "sum_to_zero")))
show("sratio cse", stancode(bf(yo ~ cse(x)), data = d, family = sratio()))
show("hurdle_cumulative equidistant", stancode(bf(yhc ~ x), data = d,
                                   family = hurdle_cumulative(threshold = "equidistant")))
show("hurdle_cumulative cs", tryCatch(stancode(bf(yhc ~ cs(x)), data = d,
                                   family = hurdle_cumulative()), error = function(e) conditionMessage(e)))
show("hurdle_cumulative thres(gr)", tryCatch(stancode(bf(yhc | thres(gr = g) ~ x), data = d,
                                   family = hurdle_cumulative()), error = function(e) conditionMessage(e)))
show("hurdle_cumulative disc~x", tryCatch(stancode(bf(yhc ~ x, disc ~ x), data = d,
                                   family = hurdle_cumulative()), error = function(e) conditionMessage(e)))
show("xbeta mixture", tryCatch(stancode(bf(yx ~ 1), data = d,
     family = mixture(xbeta, xbeta)), error = function(e) conditionMessage(e)))
show("zibb mixture", tryCatch(stancode(bf(ybb | trials(tr) ~ 1), data = d,
     family = mixture(zero_inflated_beta_binomial, beta_binomial)), error = function(e) conditionMessage(e)))
show("hc mixture", tryCatch(stancode(bf(yhc ~ 1), data = d,
     family = mixture(hurdle_cumulative, hurdle_cumulative)), error = function(e) conditionMessage(e)))
for (th in c("flexible", "equidistant", "sum_to_zero")) {
  for (f in c("cumulative", "sratio", "cratio", "acat")) {
    r <- tryCatch({brmsfamily(f, threshold = th); "ok"}, error = function(e) conditionMessage(e))
    cat(f, th, r, "\n")
  }
}
