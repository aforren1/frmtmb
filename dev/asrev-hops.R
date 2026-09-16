## BLOCKER 2 recheck: is there a THIRD hop? Read the base-R sources on
## the step()/add1()/drop1() path and list every argument stats passes
## by NAME into a generic a frmtmb method answers, then run the path.
## Usage: Rscript dev/asrev-hops.R <lib>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1],
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n\n")

cat("==== what stats passes by name into a generic ====\n")
gens <- c("nobs", "drop1", "add1", "extractAIC", "terms", "formula",
          "update", "deviance", "logLik", "model.frame", "coef",
          "vcov", "residuals", "fitted", "predict", "sigma", "weights",
          "df.residual", "family", "anova", "simulate", "confint",
          "print", "summary")
callers <- c("step", "add1.default", "drop1.default", "extractAIC.default",
             "extractAIC.lm", "sigma.default", "nobs.default",
             "AIC", "BIC", "confint.default", "profile.default",
             "update.default", "model.frame.default")
ns <- asNamespace("stats")
for (cl in callers) {
  o <- tryCatch(get(cl, envir = ns), error = function(e) NULL)
  if (!is.function(o)) { cat(sprintf("%-20s <absent>\n", cl)); next }
  src <- deparse(body(o))
  out <- character()
  for (g in gens) {
    ln <- grep(paste0("(^|[^._[:alnum:]])", g, "\\("), src, value = TRUE)
    for (l in ln) {
      nm <- regmatches(l, gregexpr("[A-Za-z._][A-Za-z._0-9]* *=",
                                   l))[[1]]
      nm <- trimws(sub("=$", "", trimws(nm)))
      if (length(nm)) out <- c(out, sprintf("%s(%s)", g,
                                            paste(nm, collapse = ",")))
    }
  }
  if (length(out)) {
    cat(sprintf("%-20s %s\n", cl,
                paste(unique(out), collapse = "  ")))
  }
}

cat("\n==== run the path ====\n")
set.seed(11)
n <- 90
dd <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                 g = factor(rep(1:9, each = 10)))
dd$y <- 0.3 + 0.5 * dd$x + rnorm(9, 0, .6)[dd$g] + rnorm(n, 0, .5)
fit <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
p <- function(lab, expr) {
  r <- tryCatch(force(expr), error = function(e) e)
  cat(sprintf("%-44s %s\n", lab,
              if (inherits(r, "error"))
                paste("ERROR:", gsub("\n", " ", conditionMessage(r)))
              else paste("ok", paste(class(r), collapse = "/"))))
}
p("nobs(fit, use.fallback = TRUE)", stats::nobs(fit, use.fallback = TRUE))
p("drop1(fit)", stats::drop1(fit))
p("drop1(fit, scale = 0, trace = 0)",
  stats::drop1(fit, scale = 0, trace = 0))
p("drop1(fit, scale = 0, trace = 0, k = 2)",
  stats::drop1(fit, scale = 0, trace = 0, k = 2))
p("extractAIC(fit)", stats::extractAIC(fit))
p("extractAIC(fit, scale = 0, k = 2)",
  stats::extractAIC(fit, scale = 0, k = 2))
p("step(fit)", utils::capture.output(stats::step(fit)))
p("step(fit, trace = 0)", utils::capture.output(stats::step(fit,
                                                            trace = 0)))
p("step(fit, trace = 0, k = log(nobs(fit)))",
  utils::capture.output(stats::step(fit, trace = 0,
                                    k = log(stats::nobs(fit)))))
p("step(fit, scope = list(lower = ~1))",
  utils::capture.output(stats::step(fit, scope = list(lower = ~ 1),
                                    trace = 0)))
p("add1(fit, ~ . + w)", stats::add1(fit, ~ . + w))
p("MuMIn::dredge-ish: terms(fit)", stats::terms(fit))
p("update(fit, . ~ . - z)", stats::update(fit, . ~ . - z))
p("AIC(fit); BIC(fit)", c(stats::AIC(fit), stats::BIC(fit)))

cat("\n==== the refusal still fires on a name stats does NOT pass ====\n")
p("drop1(fit, nosucharg = 1)", stats::drop1(fit, nosucharg = 1))
p("nobs(fit, nosucharg = 1)", stats::nobs(fit, nosucharg = 1))
p("drop1(fit, re.form = NA)", stats::drop1(fit, re.form = NA))
cat("\nDONE\n")
