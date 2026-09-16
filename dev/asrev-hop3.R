## The THIRD hop, isolated. stats::as.formula(object) is
##   formula(object, env = baseenv())
## so every caller of as.formula() on a fit reaches
## formula.frmtmb_fit(x, ...) with an `env` it now refuses.
## Usage: Rscript dev/asrev-hop3.R <lib>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1],
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n")
cat("as.formula body:\n")
print(body(stats::as.formula))

set.seed(11)
n <- 90
dd <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                 g = factor(rep(1:9, each = 10)))
dd$y <- 0.3 + 0.5 * dd$x + rnorm(9, 0, .6)[dd$g] + rnorm(n, 0, .5)
one <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
two <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
p <- function(lab, expr) {
  r <- tryCatch(force(expr), error = function(e) e)
  cat(sprintf("%-44s %s\n", lab,
              if (inherits(r, "error"))
                paste("ERROR:", gsub("\n", " ", conditionMessage(r)))
              else paste("ok", paste(class(r), collapse = "/"))))
}
cat("\n---- the one-liner ----\n")
p("as.formula(fit)", stats::as.formula(two))
p("formula(fit, env = baseenv())", stats::formula(two, env = baseenv()))
p("update.formula(fit, ~ . + w)", stats::update.formula(two, ~ . + w))
p("stats::add.scope(fit, ~ x + z + w)",
  stats::add.scope(two, stats::update.formula(two, ~ . + w)))

cat("\n---- BLOCKER 2's original fixture: ONE fixed term ----\n")
p("step(one)", utils::capture.output(stats::step(one)))
p("step(one, trace = 0)", utils::capture.output(stats::step(one,
                                                            trace = 0)))
p("nobs(one, use.fallback = TRUE)", stats::nobs(one, use.fallback = TRUE))

cat("\n---- two fixed terms ----\n")
p("step(two, trace = 0)", utils::capture.output(stats::step(two,
                                                            trace = 0)))
p("step(two, scope = list(lower = ~1), trace = 0)",
  utils::capture.output(stats::step(two, scope = list(lower = ~ 1),
                                    trace = 0)))
p("add1(two, ~ . + w)", stats::add1(two, ~ . + w))
cat("\nDONE\n")
