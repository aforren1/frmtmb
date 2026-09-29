## REVIEW claim 1, the remaining refit paths named in the findings table:
## frm_allfit(), anova(refit = TRUE), confint(method = "profile"),
## conditional_effects(band = "boot") and the autoscale pre-fit, on an
## ordinal fit whose bootstrap replicates lose the top category.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, "\n\n")
options(width = 140)

say <- function(tag, expr) {
  cat("==", tag, "\n")
  out <- tryCatch(expr, error = function(e) e)
  if (inherits(out, "condition")) {
    cat("    ERROR:", gsub("\n", " ", conditionMessage(out)), "\n\n")
  } else {
    print(out); cat("\n")
  }
}

set.seed(202)
x <- stats::rnorm(40)
cp <- cbind(stats::plogis(-0.6 - 0.5 * x), stats::plogis(0.5 - 0.5 * x),
            stats::plogis(2.6 - 0.5 * x))
dd <- data.frame(x = x, z = stats::rnorm(40),
                 y = 1L + rowSums(stats::runif(40) > cp))
cat("table(y) =", paste(table(dd$y), collapse = "/"), "\n")
fit <- suppressWarnings(frm(bf(y ~ x + z), family = cumulative(), data = dd))
cat("fitted n_tau =", length(fit$estimates[["tau_raw"]]), "\n\n")

say("frm_allfit(): threshold count of every optimizer arm",
    {
      af <- suppressWarnings(frm_allfit(fit))
      vapply(af$fits, function(f) {
        if (is.null(f)) NA_integer_ else
          length(f$estimates[["tau_raw"]])
      }, 1L)
    })

say("anova(refit = TRUE): the reduced model's threshold count",
    {
      f0 <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = dd))
      a <- suppressWarnings(anova(f0, fit, refit = TRUE))
      as.data.frame(a)[, intersect(c("npar", "Df", "Chisq", "logLik"),
                                   names(as.data.frame(a)))]
    })

say("confint(method = 'profile') on the top threshold",
    suppressWarnings(confint(fit, method = "profile",
                             parm = "tau_raw_3")))

say("conditional_effects(band = 'boot')",
    {
      ce <- suppressWarnings(conditional_effects(fit, effects = "x",
                                                band = "boot", ndraws = 20,
                                                seed = 31))
      d <- as.data.frame(ce[[1L]])
      c(rows = nrow(d), na_lower = sum(is.na(d[["lower__"]])),
        na_upper = sum(is.na(d[["upper__"]])))
    })

say("the autoscale pre-fit keeps the count",
    {
      fa <- suppressWarnings(frm(bf(y ~ x + z), family = cumulative(),
                                data = dd,
                                control = frmtmb_control(autoscale = TRUE)))
      c(n_tau = length(fa$estimates[["tau_raw"]]))
    })
cat("DONE rev-15\n")
