# Smoke test of the lane build: disc and threshold structures on the
# ordinal families. Output: dev/ordinal-log-smoke1.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(20260930)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
try1 <- function(label, expr) {
  cat("\n==", label, "==\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    NULL
  })
  invisible(r)
}
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  for (th in c("flexible", "equidistant", "sum_to_zero")) {
    f <- try1(paste(fam, th, "disc ~ 0 + z"), {
      frm(bf(y ~ x, disc ~ 0 + z),
          family = get(fam)(threshold = th), data = d)
    })
    if (is.null(f)) next
    cat("logLik", format(as.numeric(logLik(f)), digits = 10), "\n")
    print(variables(f))
    try1("fitted dims", print(dim(fitted(f))))
    try1("fitted linear", print(head(fitted(f, scale = "linear"), 2)))
    try1("predict", print(head(predict(f), 2)))
    try1("simulate", print(table(simulate(f, nsim = 1, seed = 1)[[1]])))
  }
}
f <- try1("summary cumulative equidistant", {
  frm(bf(y ~ x), family = cumulative(threshold = "equidistant"), data = d)
})
print(summary(f))
print(fixef(f))
print(default_prior(bf(y ~ x), family = cumulative(threshold = "equidistant"),
                    data = d))
print(hypothesis(f, "delta > 0.5", class = NULL))
try1("hurdle equidistant", {
  fh <- frm(bf(yh ~ x), family = hurdle_cumulative(threshold = "equidistant"),
            data = d)
  print(variables(fh)); print(summary(fh))
})
try1("hurdle sum_to_zero", {
  fh <- frm(bf(yh ~ x), family = hurdle_cumulative(threshold = "sum_to_zero"),
            data = d)
  print(variables(fh)); print(fitted(fh)[1:2, , ])
})
try1("grouped equidistant", {
  fg <- frm(bf(y | thres(gr = g) ~ x),
            family = sratio(threshold = "equidistant"), data = d)
  print(variables(fg)); print(summary(fg))
})
try1("grouped stz", {
  fg <- frm(bf(y | thres(gr = g) ~ x),
            family = cumulative(threshold = "sum_to_zero"), data = d)
  print(variables(fg))
})
try1("disc intercept warns", {
  frm(bf(y ~ x, disc ~ z), family = acat(), data = d)
})
try1("delta prior", {
  fp <- frm(bf(y ~ x), family = cumulative(threshold = "equidistant"),
            data = d, prior = set_prior("normal(0, 1)", class = "delta"))
  print(fp$prior)
  print(variables(fp))
})
try1("stz intercept prior refused", {
  frm(bf(y ~ x), family = cumulative(threshold = "sum_to_zero"),
      data = d, prior = set_prior("normal(0, 1)", class = "Intercept"))
})
