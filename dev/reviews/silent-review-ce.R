# Review of item 720: no false refusal on an evaluable vreal() / rate().
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
show <- function(lab, expr) {
  r <- tryCatch(suppressWarnings(expr),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r)) cat("==", lab, substr(r, 1, 200), "\n") else {
    df <- r[[1L]]
    cat("==", lab, " rows", nrow(df), " estimate range",
        paste(format(range(df$estimate__), digits = 5), collapse = " to "),
        " NA", sum(is.na(df$estimate__)), "\n")
  }
}
scaled <- custom_family(
  "rv_scaled", dpars = c("mu", "sigma"),
  links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms) {
    RTMB::dnorm(y, dpars[["mu"]] * aterms[["vreal1"]], dpars[["sigma"]],
                log = TRUE)
  },
  post = list(mean_fn = function(dpars, aterms) {
    dpars[["mu"]] * aterms[["vreal1"]]
  }),
  sim = function(dpars, aterms, n) {
    stats::rnorm(n, dpars[["mu"]] * aterms[["vreal1"]], dpars[["sigma"]])
  })
set.seed(3)
n <- 120
d <- data.frame(x = rnorm(n), s = rep_len(c(1, 3), n),
                t = runif(n, 1, 5))
d$y <- (1 + 0.5 * d$x) * d$s + rnorm(n, 0, 0.3)
d$cnt <- rpois(n, d$t * exp(0.2 + 0.3 * d$x))
fv <- frm(bf(y | vreal(s) ~ x), family = scaled, data = d)
fr <- frm(bf(cnt | rate(t) ~ x), family = poisson(), data = d)
show("core vreal mean", conditional_effects(fv, effects = "x"))
show("core vreal predict", conditional_effects(fv, effects = "x",
                                               method = "predict"))
show("core rate mean", conditional_effects(fr, effects = "x"))
show("core rate predict", conditional_effects(fr, effects = "x",
                                              method = "predict"))
sv <- tryCatch(frm_sample(fv, chains = 1, iter = 300, warmup = 150,
                          seed = 1, refresh = 0),
               error = function(e) conditionMessage(e))
if (is.character(sv)) cat("sample vreal ERROR", sv, "\n") else {
  show("sample vreal predict", conditional_effects(sv, effects = "x",
                                                   method = "predict"))
  show("sample vreal mean", conditional_effects(sv, effects = "x"))
}
sr <- tryCatch(frm_sample(fr, chains = 1, iter = 300, warmup = 150,
                          seed = 1, refresh = 0),
               error = function(e) conditionMessage(e))
if (is.character(sr)) cat("sample rate ERROR", sr, "\n") else {
  show("sample rate predict", conditional_effects(sr, effects = "x",
                                                  method = "predict"))
}
