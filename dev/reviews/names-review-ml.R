# Review of lane names: ML-fit readers on both arms.
# Usage: Rscript dev/reviews/names-review-ml.R <base|lane>
# Output: dev/reviews/names-review-ml-<arm>.rds and stdout.
a <- commandArgs(TRUE)
arm <- a[1]
lib <- c(if (arm == "lane") "C:/Users/adf44/source/r/wt-names-lib",
         "C:/Users/adf44/source/r/rellib-r20",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
res <- list()
try_out <- function(label, expr) {
  v <- tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                   class = "err"))
  cat("\n==== ", label, "\n", sep = "")
  if (inherits(v, "err")) cat("ERROR:", v, "\n") else print(v)
  res[[label]] <<- v
  invisible(v)
}
set.seed(11)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(rep(sprintf("g%02d", 1:12), each = 20)),
                h = factor(rep(1:8, length.out = n)),
                f = factor(rep(c("fa", "fb", "fc"), length.out = n)))
ug <- rnorm(12, 0, 0.6)[d$g]
us <- rnorm(12, 0, 0.4)[d$g]
d$y <- 1 + 0.5 * d$x + ug + us * d$x + rnorm(n, 0, exp(0.2 * d$z + 0.3 * us))
d$y2 <- 0.5 * d$y + rnorm(n)
q <- function(e) suppressWarnings(suppressMessages(e))

f_two <- q(frm(y ~ x + (1 | g) + (0 + x | g), data = d))
f_dist <- q(frm(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), data = d))
f_mv <- q(frm(bf(y ~ x + (1 | g)) + bf(y2 ~ z + (1 | g)) + set_rescor(FALSE),
              data = d))
f_three <- q(frm(bf(y ~ x + (1 + x | g) + (1 | h), sigma ~ (1 | g)), data = d))

for (nm in c("f_two", "f_dist", "f_mv", "f_three")) {
  fit <- get(nm)
  try_out(paste(nm, "variables"), variables(fit))
  try_out(paste(nm, "ranef names"), names(ranef(fit)))
  re <- ranef(fit, condVar = TRUE)
  try_out(paste(nm, "ranef cols"), lapply(unclass(re), colnames))
  try_out(paste(nm, "ranef long"), {
    df <- as.data.frame(re)
    df[order(df$grp, df$term, df$level), ]
  })
  try_out(paste(nm, "coef"), lapply(coef(fit), function(x) {
    if (is.data.frame(x)) head(x, 3) else lapply(x, head, 3)
  }))
  try_out(paste(nm, "VarCorr names"), names(VarCorr(fit)))
  try_out(paste(nm, "confint_varcorr"), confint_varcorr(fit)[, 1:4])
  try_out(paste(nm, "insight find_random"), insight::find_random(fit))
  try_out(paste(nm, "insight get_parameters"), insight::get_parameters(fit))
  try_out(paste(nm, "insight get_variance"),
          q(insight::get_variance(fit)))
  try_out(paste(nm, "insight get_random"),
          head(q(insight::get_random(fit))))
  try_out(paste(nm, "insight find_parameters"),
          q(insight::find_parameters(fit)))
  try_out(paste(nm, "marginaleffects avg_predictions"),
          q(as.data.frame(marginaleffects::avg_predictions(fit)))[, 1:3])
  try_out(paste(nm, "emmeans"),
          q(as.data.frame(emmeans::emmeans(fit, ~ x, at = list(x = 0)))))
  try_out(paste(nm, "lme4-style $g[1:2,]"), re$g[1:2, , drop = FALSE])
}
try_out("f_dist [[sigma: 1 | g]]", ranef(f_dist)[["sigma: 1 | g"]][1:3, , drop = FALSE])
try_out("f_dist [[1 | g]]", ranef(f_dist)[["1 | g"]][1:3, , drop = FALSE])
try_out("f_two [[0 + x | g]]", ranef(f_two)[["0 + x | g"]][1:3, , drop = FALSE])
try_out("f_two print", capture.output(print(ranef(f_two)))[1:4])
try_out("f_dist ranef groups=", names(ranef(f_dist, groups = "g")))
try_out("f_dist ranef pars=", lapply(unclass(ranef(f_dist, pars = "sigma_Intercept")), colnames))
try_out("f_three ranef pars= Intercept",
        lapply(unclass(ranef(f_three, pars = "Intercept")), colnames))

# no group terms
f_s <- q(frm(y ~ s(x), data = d))
f_none <- q(frm(y ~ x, data = d))
try_out("smooth ranef", ranef(f_s))
try_out("none ranef", ranef(f_none))
try_out("none coef", coef(f_none))
try_out("smooth coef", coef(f_s))
try_out("smooth insight find_random", insight::find_random(f_s))
try_out("smooth insight get_variance", q(insight::get_variance(f_s)))

# || display label
f_db <- q(frm(y ~ x + (1 + f || g), data = d))
try_out("dblbar confint_varcorr", confint_varcorr(f_db)[, 1:4])
try_out("dblbar VarCorr", names(VarCorr(f_db)))
f_gr <- q(frm(y ~ x + (1 + x | gr(g, cor = FALSE)), data = d))
try_out("gr cor=FALSE confint_varcorr", confint_varcorr(f_gr)[, 1:4])

# summary Formula and Data
try_out("summary dist head", capture.output(summary(f_dist))[1:6])
try_out("summary mv head", capture.output(summary(f_mv))[1:6])
f_nl <- q(frm(bf(y ~ a + b * x, a ~ 1 + (1 | g), b ~ 1, nl = TRUE), data = d))
try_out("summary nl head", capture.output(summary(f_nl))[1:7])
try_out("variables nl", variables(f_nl))
f_fix <- q(frm(bf(y ~ x, nu = 4), family = student(), data = d))
try_out("summary pfix head", capture.output(summary(f_fix))[1:6])
f_long <- q(frm(bf(y ~ x + z + I(x^2) + I(z^2) + I(x * z) + I(x^3) + I(z^3) +
                   I(x^2 * z) + I(x * z^2) + (1 | g),
                 sigma ~ z + x + I(x^2) + I(z^2) + I(x * z) + I(x^3) +
                   I(z^3) + I(x^2 * z)), data = d))
try_out("summary long formula", capture.output(summary(f_long))[1:8])
try_out("print long formula", capture.output(print(f_long))[1:8])

# predict(summary = FALSE) dimnames
try_out("predict summary=FALSE dimnames",
        dimnames(predict(f_dist, summary = FALSE, ndraws = 5)))
saveRDS(res, sprintf("dev/reviews/names-review-ml-%s.rds", arm))
cat("\nDONE\n")
