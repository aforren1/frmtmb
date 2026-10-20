# Lane sampleprior, fix round 1, item 2: what brms 2.23.0 does with a
# bounded density on a parameter that has no bound of its own, class by
# class (data seed 20261021). brms's check_prior_content() is what
# speaks; the Stan declaration shows whether brms bounds the parameter.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
set.seed(20261021)
d <- data.frame(x = stats::rnorm(40), g = factor(rep(1:5, 8)))
d$y <- 1 + 2.4 * d$x + stats::rnorm(5)[d$g] + stats::rnorm(40)
cases <- list(
  "uniform on b" = brms::set_prior("uniform(2.2, 2.6)", class = "b",
                                   coef = "x"),
  "uniform on b, lb ub" = brms::set_prior("uniform(2.2, 2.6)", class = "b",
                                          lb = 2.2, ub = 2.6),
  "gamma on b" = brms::set_prior("gamma(2, 1)", class = "b", coef = "x"),
  "beta on Intercept" = brms::set_prior("beta(2, 2)", class = "Intercept"),
  "uniform on sd" = brms::set_prior("uniform(0, 10)", class = "sd"),
  "uniform on sigma" = brms::set_prior("uniform(0, 10)", class = "sigma"),
  "exponential on sd" = brms::set_prior("exponential(1)", class = "sd"),
  "lognormal on sigma" = brms::set_prior("lognormal(0, 1)", class = "sigma"))
for (nm in names(cases)) {
  # brms takes lb and ub on a whole class only ("coef may not be
  # specified when using boundaries")
  w <- character(0)
  code <- withCallingHandlers(
    brms::stancode(y ~ x + (1 | g), data = d, prior = cases[[nm]]),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  s <- strsplit(as.character(code), "\n")[[1]]
  decl <- grep("^ *(vector|real).*(  b;| b;|Intercept;|sigma;|sd_1;)", s,
               value = TRUE)
  cat(sprintf("== %s: %d warning(s)\n", nm, length(w)))
  for (x in w) cat("  WARNING:", gsub("\n", " | ", x), "\n")
  cat("  ", trimws(decl), sep = "\n  ")
  cat("\n")
}
