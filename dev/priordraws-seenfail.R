# The gaps the new tests pin, constructed with arguments rellib-r14
# already takes, so the base arm fails on behavior and not on a missing
# argument. Usage: Rscript dev/priordraws-seenfail.R <lane|base>
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (identical(arm, "lane")) "C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "lib:", find.package("frmtmb.sample"),
    format(packageVersion("frmtmb.sample")), "\n")
show <- function(lab, expr) {
  cat("\n--", lab, "\n")
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n"); invisible(NULL)
    }),
    warning = function(w) {
      cat("WARNING:", substr(conditionMessage(w), 1, 160), "\n")
      invokeRestart("muffleWarning")
    })
  if (!is.null(r)) print(r)
}
set.seed(21)
n <- 150
dd <- data.frame(x = rnorm(n), w = rnorm(n), g = factor(rep(1:10, 15)))
dd$y <- rnorm(n, 0.5 + 0.6 * dd$w + rnorm(10, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + w + (1 | g)), family = gaussian(), data = dd)
ds <- suppressWarnings(suppressMessages(frm_sample(
  fit, chains = 2, iter = 1200, warmup = 400, refresh = 0, seed = 7,
  prior = set_prior("normal(0, 1)", class = "b"))))
show("A dim(as_draws_array(inc_warmup = TRUE)), brms: 1200 x 2",
     dim(as_draws_array(ds, inc_warmup = TRUE)))
show("B dim(as.mcmc(inc_warmup = TRUE)[[1]])[1], brms: 1200",
     dim(as.mcmc(ds, inc_warmup = TRUE)[[1]])[1])
show("C Evid.Ratio of x = 0 without prior draws, brms: NA, no warning",
     hypothesis(ds, "x = 0")$hypothesis$Evid.Ratio)
show("D update(ds, sample_prior = 'yes') prior columns, brms: b, Intercept, sigma, sd_g",
     grep("^prior_", variables(update(ds, sample_prior = "yes")),
          value = TRUE))
show("E brms::prior_draws(ds) dim, brms on a sample_prior = 'yes' fit: 1600 rows",
     dim(brms::prior_draws(ds)))
show("F Evid.Ratio of sigma = 1 on draws with prior draws, brms: a number",
     {
       d2 <- tryCatch(suppressWarnings(suppressMessages(
         update(ds, sample_prior = "yes"))), error = function(e) ds)
       hypothesis(d2, "sigma = 1", class = NULL)$hypothesis$Evid.Ratio
     })
