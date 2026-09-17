# Reviewer (brmsport): when did the prior table start validating the
# response against the family? One commit per process, loaded with
# pkgload from a `git archive` extract in the session scratchpad.
#   Rscript dev/brmsport-rev-regress.R <srcdir> <label>
# Data seed 91.
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
a <- commandArgs(trailingOnly = TRUE)
suppressMessages(pkgload::load_all(a[1], export_all = FALSE, quiet = TRUE,
                                   helpers = FALSE, attach_testthat = FALSE))
show <- function(label, expr) {
  out <- tryCatch(paste(format(expr), collapse = " "),
                  error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-10s %-44s %s\n", a[2], label, substr(out, 1, 120)))
}
set.seed(91)
dat <- data.frame(y = rnorm(10), x = rnorm(10), z = rnorm(10),
                  g = rep(1:2, 5))
db <- dat; db$y <- plogis(db$y)
dn <- data.frame(y = rbinom(10, 5, 0.5), x = rnorm(10))
has_dp <- exists("default_prior", mode = "function")
tab <- function(f, d, fam = NULL) {
  if (has_dp) nrow(default_prior(f, data = d, family = fam))
  else nrow(get_prior(f, data = d, family = fam))
}
cat(a[2], "version", as.character(utils::packageVersion("frmtmb")),
    "default_prior exists:", has_dp, "\n")
show("get_prior Beta, y = rnorm (brms: 7 rows)",
     nrow(get_prior(bf(y ~ 1, phi ~ z + (1 | g)) + Beta(), data = dat)))
show("table Beta, y = rnorm, family in bf()",
     tab(bf(y ~ 1, phi ~ z + (1 | g), family = Beta()), dat))
show("table Beta, y in (0, 1) control",
     tab(bf(y ~ 1, phi ~ z + (1 | g), family = Beta()), db))
show("get_prior binomial, no trials()",
     nrow(get_prior(y ~ x, data = dn, family = binomial())))
show("get_prior poisson, y = rnorm",
     nrow(get_prior(y ~ x, data = dat, family = poisson())))
