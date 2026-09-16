## A FOURTH hop? The class-agnostic survivors of dev/asrev-hopscan.R,
## run rather than counted, change against base.
## Usage: Rscript dev/asrev-hop4.R <lib>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1],
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(11)
n <- 90
dd <- data.frame(x = rnorm(n), z = rnorm(n),
                 g = factor(rep(1:9, each = 10)))
dd$y <- 0.3 + 0.5 * dd$x + rnorm(9, 0, .6)[dd$g] + rnorm(n, 0, .5)
fit <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
p <- function(lab, expr) {
  r <- tryCatch(force(expr), error = function(e) e)
  cat(sprintf("%-42s %s\n", lab,
              if (inherits(r, "error"))
                paste("ERROR:", gsub("\n", " ", conditionMessage(r)))
              else paste("ok", paste(class(r), collapse = "/"))))
}
p("as.formula(fit)", stats::as.formula(fit))
p("terms(fit)", stats::terms(fit))
p("terms(fit, data = dd)", stats::terms(fit, data = dd))
p("model.matrix(fit)", stats::model.matrix(fit))
p("model.frame(fit)", stats::model.frame(fit))
p("model.frame(fit, data = dd)", stats::model.frame(fit, data = dd))
p("influence.measures(fit)", stats::influence.measures(fit))
p("profile(fit)", stats::profile(fit))
p("confint(fit, trace = FALSE)", stats::confint(fit, trace = FALSE))
p("case.names(fit)", stats::case.names(fit))
p("variable.names(fit)", stats::variable.names(fit))
p("formula(fit)", stats::formula(fit))
p("getCall(fit)", stats::getCall(fit))
p("effects(fit)", stats::effects(fit))
p("proj(fit)", stats::proj(fit))
p("dummy.coef(fit)", stats::dummy.coef(fit))
p("relevel-ish: update(fit, data = dd)", stats::update(fit, data = dd))
p("simulate(fit, nsim = 2, seed = 1)",
  stats::simulate(fit, nsim = 2, seed = 1))
p("residuals(fit, na.rm = TRUE)", stats::residuals(fit, na.rm = TRUE))
p("df.residual(fit, na.rm = TRUE)", stats::df.residual(fit, na.rm = TRUE))
p("coef(fit, digits = 3)", stats::coef(fit, digits = 3))
cat("\nDONE\n")
