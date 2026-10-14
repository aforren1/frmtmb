# Smoke test: every path for the three new families on one data set each.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
try1 <- function(label, expr) {
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  WARN [", label, "]:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) {
    cat("  ERROR [", label, "]:", conditionMessage(e), "\n")
    NULL
  })
  invisible(r)
}
set.seed(20261013)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(20, n / 20))
re <- rnorm(20, 0, 0.3)[d$g]
nu <- exp(log(3) + 0.2 * d$z) + 1
mu_f <- exp(1 + 0.3 * d$x + re)
d$yf <- mu_f / gamma(1 - 1 / nu) * rexp(n)^(-1 / nu)
xi <- 0.15
d$yg <- 1 + 0.5 * d$x + re + exp(-0.3 + 0.2 * d$z) * (rexp(n)^(-xi) - 1) / xi
q <- plogis(1 + 0.4 * d$x + re)
k <- exp(0.2 + 0.2 * d$z)
d$yd <- pmax(ceiling((log1p(-runif(n)) / log(q))^(1 / k) - 1), 0)

specs <- list(
  frechet = list(f = bf(yf ~ x + (1 | g), nu ~ z), fam = frechet()),
  gev = list(f = bf(yg ~ x + (1 | g), sigma ~ z), fam = gen_extreme_value()),
  dw = list(f = bf(yd ~ x + (1 | g), shape ~ z), fam = discrete_weibull())
)
for (nm in names(specs)) {
  cat("\n=====", nm, "\n")
  s <- specs[[nm]]
  fit <- try1("frm", frm(s$f, family = s$fam, data = d))
  if (is.null(fit)) next
  print(summary(fit))
  try1("fitted", print(head(fitted(fit))))
  try1("predict", print(head(predict(fit))))
  try1("frm_linpred response", print(head(frm_linpred(fit, type = "response"))))
  try1("residuals", print(head(residuals(fit))))
  try1("residuals pearson", print(head(residuals(fit, type = "pearson"))))
  try1("dharma",
       print(class(dharma_residuals(fit, n = 50))))
  try1("residuals osa", print(head(residuals(fit, type = "osa"))))
  try1("simulate", print(dim(simulate(fit, nsim = 3, seed = 1))))
  try1("log_lik", print(dim(log_lik(fit))))
  try1("posterior_predict", print(dim(posterior_predict(fit, ndraws = 5))))
  try1("posterior_epred", print(dim(posterior_epred(fit, ndraws = 5))))
  try1("conditional_effects", print(class(conditional_effects(fit))))
  try1("pp_check", print(class(pp_check(fit, ndraws = 5))))
  try1("default_prior", print(default_prior(s$f, data = d, family = s$fam)))
  try1("logLik", print(logLik(fit)))
}
cat("\n===== cens/trunc\n")
d$cen <- ifelse(d$yf > quantile(d$yf, 0.8), 1, 0)
d$yfc <- pmin(d$yf, quantile(d$yf, 0.8))
try1("frechet cens", print(fixef(frm(yfc | cens(cen) ~ x, family = frechet(),
                                     data = d))))
try1("frechet trunc", print(fixef(frm(yf | trunc(lb = 0.5) ~ x,
                                      family = frechet(),
                                      data = d[d$yf > 0.5, ]))))
d$cg <- ifelse(d$yg > quantile(d$yg, 0.85), 1, 0)
d$ygc <- pmin(d$yg, quantile(d$yg, 0.85))
try1("gev cens", print(fixef(frm(ygc | cens(cg) ~ x,
                                 family = gen_extreme_value(), data = d))))
d$cd <- ifelse(d$yd >= 6, 1, 0)
d$ydc <- pmin(d$yd, 6)
try1("dw cens", print(fixef(frm(ydc | cens(cd) ~ x,
                                family = discrete_weibull(), data = d))))
dt <- d[d$yd >= 1, ]
fdt <- try1("dw trunc", frm(yd | trunc(lb = 1) ~ x,
                            family = discrete_weibull(), data = dt))
if (!is.null(fdt)) {
  print(fixef(fdt))
  try1("dw trunc fitted", print(head(fitted(fdt))))
}
cat("\n===== Gumbel via xi = 0\n")
d$yu <- 1 + 0.5 * d$x - log(rexp(n))
try1("gumbel", print(summary(frm(bf(yu ~ x, xi = 0),
                                 family = gen_extreme_value(), data = d))))
cat("\n===== mixtures\n")
d$ym <- ifelse(runif(n) < 0.5, rpois(n, 1), rpois(n, 12))
try1("dw mixture", print(summary(frm(bf(ym ~ 1),
  family = mixture(discrete_weibull(), discrete_weibull()), data = d))))
cat("\n===== string families\n")
try1("frechet string", print(class(frm(yf ~ x, family = "frechet", data = d))))
try1("wiener string", frm(yf ~ x, family = "wiener", data = d))
