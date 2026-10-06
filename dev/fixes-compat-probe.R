# Lane fixes, item 4d: the cells the compat rows for disc and the two
# threshold structures will claim, measured.
#   Rscript dev/fixes-compat-probe.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
set.seed(62)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)),
                h = factor(rep(1:10, 20)))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yc <- rnorm(n)
try_ <- function(lab, expr) {
  r <- tryCatch({
    v <- suppressWarnings(suppressMessages(expr))
    paste("OK", paste(class(v)[1], collapse = ""))
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 160)))
  cat(sprintf("%-40s %s\n", lab, r))
}
try_("gaussian disc ~ x", frm(bf(yc ~ x, disc ~ z), data = d))
try_("gaussian(threshold =)", gaussian(threshold = "equidistant"))
try_("poisson(threshold =)", poisson(threshold = "equidistant"))
nd <- d[1:5, ]
cases <- list(
  disc = list(bf(y ~ x, disc ~ 0 + z), cumulative(),
              bf(y ~ x + (1 | h), disc ~ 0 + z),
              bf(y | thres(gr = g) ~ x, disc ~ 0 + z)),
  equidistant = list(bf(y ~ x), cumulative(threshold = "equidistant"),
                     bf(y ~ x + (1 | h)), bf(y | thres(gr = g) ~ x)),
  sum_to_zero = list(bf(y ~ x), sratio(threshold = "sum_to_zero"),
                     bf(y ~ x + (1 | h)), bf(y | thres(gr = g) ~ x)))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- frm(cs[[1]], family = cs[[2]], data = d)
  cat("==", nm, "\n")
  try_("fitted", fitted(fit))
  try_("predict(newdata)", predict(fit, newdata = nd))
  try_("simulate", simulate(fit, nsim = 2, seed = 1))
  try_("residuals osa", residuals(fit, type = "osa"))
  try_("emmeans", emmeans::emmeans(fit, ~ x))
  try_("REML", frm(cs[[1]], family = cs[[2]], data = d, REML = TRUE))
  try_("random intercept + quadrature",
       frm(cs[[3]], family = cs[[2]], data = d, quadrature = TRUE))
  try_("thres(gr = g)", frm(cs[[4]], family = cs[[2]], data = d))
  try_("mixture", frm(cs[[1]], family = mixture(cs[[2]], cs[[2]]),
                      data = d))
  try_("mvbf with gaussian",
       frm(mvbf(cs[[1]], bf(yc ~ x)), family = list(cs[[2]], gaussian()),
           data = d))
  try_("confint profile", confint(fit, parm = "x", method = "profile"))
}
