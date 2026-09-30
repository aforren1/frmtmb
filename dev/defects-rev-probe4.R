# Reviewer of lane defects: mixture(order =) identity, hsgp lscale scale,
# brms on a continuous response to cumulative().
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
a <- mixture(gaussian(), gaussian()); b <- mixture(gaussian(), gaussian())
cat("control: two NULL calls identical:", identical(a, b), "\n")
cat("all.equal NULL vs none:", isTRUE(all.equal(a, mixture(gaussian(), gaussian(), order = "none"))), "\n")
cat("identical ignoring environments:", identical(a, mixture(gaussian(), gaussian(), order = "none"), ignore.environment = TRUE), "\n")
set.seed(141)
n <- 150
d <- data.frame(x = runif(n, 0, 10))
d$y <- sin(d$x) + rnorm(n, 0, 0.3)
for (fo in list(y ~ gp(x), y ~ gp(x, k = 20), y ~ gp(x, k = 40))) {
  f <- frm(fo, data = d)
  s <- summary(f)$gp
  cat(deparse1(fo), " logLik", format(as.numeric(logLik(f))), "\n")
  print(s[, 1:2])
  bk <- Filter(function(b) b$covstruct %in% c("gp", "hsgp"), f$frame$re_blocks)[[1]]
  cat("  covstruct", bk$covstruct, " gp_dmax", format(bk$gp_dmax %||% NA), " has aux_D2", !is.null(bk$aux_D2), "\n")
  print(tryCatch(confint_varcorr(f), error = function(e) conditionMessage(e)))
}
suppressMessages(library(brms))
dd <- data.frame(x = rnorm(30), yb = rnorm(30))
print(tryCatch(brms::default_prior(yb ~ x, data = dd, family = brms::cumulative()),
               error = function(e) conditionMessage(e)))
print(tryCatch(brms::validate_prior(brms::set_prior("normal(0,1)", class = "b"), yb ~ x, data = transform(dd, yb = yb + 5), family = brms::Beta()),
               error = function(e) conditionMessage(e)))
