## REVIEW 2.5e: hunt a call a GENERIC's own contract permits that the
## new refusal now rejects. Compared arm by arm against 0.57.0.
## Usage: Rscript dev/asrev-overrefuse.R <lib>
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n")

set.seed(11)
n <- 90
dd <- data.frame(x = rnorm(n), g = factor(rep(1:9, each = 10)))
dd$y <- 0.3 + 0.5 * dd$x + rnorm(9, 0, .6)[dd$g] + rnorm(n, 0, .5)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

p <- function(lab, expr) {
  r <- tryCatch(force(expr), error = function(e) e)
  cat(sprintf("%-52s %s\n", lab,
              if (inherits(r, "error"))
                paste("ERROR:", gsub("\n", " ", conditionMessage(r)))
              else paste("ok", paste(class(r), collapse = "/"))))
}

cat("\n-- stats generics whose documented contract carries extra args --\n")
p("model.frame(fit, data = dd)", stats::model.frame(fit, data = dd))
p("model.frame(fit)", stats::model.frame(fit))
p("print(fit, digits = 3)", utils::capture.output(print(fit, digits = 3)))
p("summary(fit, digits = 3)", summary(fit, digits = 3))
p("print(summary(fit), digits = 3)",
  utils::capture.output(print(summary(fit), digits = 3)))
p("logLik(fit, REML = TRUE)", stats::logLik(fit, REML = TRUE))
p("nobs(fit, use.fallback = TRUE)", stats::nobs(fit, use.fallback = TRUE))
p("family(fit, ...)", stats::family(fit))
p("formula(fit, env = globalenv())", stats::formula(fit, env = globalenv()))
p("coef(fit, complete = TRUE)", stats::coef(fit, complete = TRUE))
p("vcov(fit, complete = TRUE)", stats::vcov(fit, complete = TRUE))
p("confint(fit, parm = 'x', level = 0.9)",
  stats::confint(fit, parm = "x", level = 0.9))
p("simulate(fit, nsim = 2, seed = 1)",
  stats::simulate(fit, nsim = 2, seed = 1))
p("residuals(fit, type = 'pearson')",
  stats::residuals(fit, type = "pearson"))
p("terms(fit)", stats::terms(fit))
p("AIC(fit)", stats::AIC(fit))
p("BIC(fit)", stats::BIC(fit))
p("sigma(fit)", stats::sigma(fit))
p("df.residual(fit)", stats::df.residual(fit))
p("update(fit, . ~ x)", stats::update(fit, . ~ x))
p("weights(fit)", stats::weights(fit))
p("str(fit) (autoprint path)", { utils::capture.output(utils::str(fit)); TRUE })
p("format(fit)", tryCatch(format(fit), error = function(e) stop(e)))

cat("\n-- R's own machinery --\n")
p("do.call(fitted, list(fit))", do.call(stats::fitted, list(fit)))
p("do.call('fitted', list(fit, ndraws = 2))",
  do.call("fitted", list(fit, ndraws = 2)))
p("lapply(list(fit), fitted)[[1]]", lapply(list(fit), stats::fitted)[[1]])
p("Map(fitted, list(fit))[[1]]", Map(stats::fitted, list(fit))[[1]])

cat("\n-- the message frm_check_dots produces from a few shapes --\n")
for (e in list(quote(fitted(fit, ndraws = 2)),
               quote(fitted(fit, re_frmula = NA)),
               quote(nobs(fit, 7)),
               quote(predict(fit, re.form = NA)),
               quote(predict(fit, allow.new.levels = TRUE)),
               quote(coef(fit, robust = TRUE)),
               quote(print(fit, digits = 3)),
               quote(model.frame(fit, data = dd)))) {
  m <- tryCatch({ eval(e); "NO ERROR" }, error = conditionMessage)
  cat(sprintf("  %-40s %s\n", deparse1(e), gsub("\n", " ", m)))
}
cat("\nDONE\n")
