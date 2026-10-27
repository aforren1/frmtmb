# brms 2.23.0's default prior, Stan declaration and prior statement for
# each autocorrelation class. No compile: default_prior() and stancode().
# Log: dev/acprior-out/brms-defaults.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 40
d <- data.frame(y = rnorm(n), x = rnorm(n), t = rep(1:10, 4),
                g = rep(1:4, each = 10))
cases <- list(
  ar1 = y ~ x + ar(time = t, gr = g, p = 1),
  ar2 = y ~ x + ar(time = t, gr = g, p = 2),
  ma1 = y ~ x + ma(time = t, gr = g),
  arma = y ~ x + arma(time = t, gr = g),
  ar1cov = y ~ x + ar(time = t, gr = g, p = 1, cov = TRUE),
  ma1cov = y ~ x + ma(time = t, gr = g, cov = TRUE),
  armacov = y ~ x + arma(time = t, gr = g, cov = TRUE),
  cosy = y ~ x + cosy(time = t, gr = g),
  unstr = y ~ x + unstr(time = t, gr = g)
)
for (nm in names(cases)) {
  cat("\n=====", nm, ":", deparse(cases[[nm]]), "\n")
  dp <- default_prior(cases[[nm]], data = d)
  print(as.data.frame(dp)[, c("prior", "class", "coef", "group", "resp",
                              "dpar", "lb", "ub", "source")])
  sc <- stancode(cases[[nm]], data = d)
  ln <- strsplit(sc, "\n")[[1]]
  hit <- grep("ar;|ma;|cosy;|cortime|Lcortime| ar\\[| ma\\[|target \\+=.*(ar|ma|cosy|cortime)|lprior \\+=.*(ar|ma|cosy|cortime)",
              ln, value = TRUE)
  cat(paste0("  | ", hit), sep = "\n")
}

# a user prior: where brms puts it and whether a bound comes with it
cat("\n===== user priors\n")
up <- list(
  ar1_normal = list(cases$ar1, set_prior("normal(0.5, 0.2)", class = "ar")),
  ar1cov_normal = list(cases$ar1cov, set_prior("normal(0.5, 0.2)",
                                               class = "ar")),
  ar1_lbub = list(cases$ar1, set_prior("normal(0.5, 0.2)", class = "ar",
                                       lb = 0, ub = 1)),
  cosy_beta = list(cases$cosy, set_prior("beta(2, 2)", class = "cosy")),
  ma1_normal = list(cases$ma1, set_prior("normal(0, 0.5)", class = "ma"))
)
for (nm in names(up)) {
  cat("\n---", nm, "\n")
  sc <- stancode(up[[nm]][[1]], data = d, prior = up[[nm]][[2]])
  ln <- strsplit(sc, "\n")[[1]]
  hit <- grep("ar;|ma;|cosy;|lprior \\+=.*(ar|ma|cosy)", ln, value = TRUE)
  cat(paste0("  | ", hit), sep = "\n")
}
# a prior on theta or on an internal name has no brms class
cat("\n--- class theta on ar:\n")
r <- tryCatch(validate_prior(set_prior("normal(0,1)", class = "theta"),
                             cases$ar1, data = d),
              error = function(e) paste("ERROR:", conditionMessage(e)))
print(r)
