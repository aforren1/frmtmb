.libPaths(c("C:/Users/adf44/source/r/rellib-r9", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
set.seed(1); d <- data.frame(x = rnorm(60)); d$y <- cut(d$x + rlogis(60), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
s <- strsplit(as.character(brms::stancode(y ~ x, data = d, family = brms::cumulative(threshold = "equidistant"))), "\n")[[1]]
cat(grep("delta", s, value = TRUE), sep = "\n")
w <- character(0)
invisible(withCallingHandlers(brms::validate_prior(brms::set_prior("lognormal(0, 1)", class = "delta"), y ~ x, data = d,
  family = brms::cumulative(threshold = "equidistant")), warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }))
cat("brms warnings on lognormal(0,1) class delta:", length(w), "\n")
