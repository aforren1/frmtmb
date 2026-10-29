# Lane meme: what frmtmb does with me() today. Run against rellib-r15
# (the "before" arm) or against the lane library, by the first argument:
#   Rscript dev/meme-before.R base > dev/meme-out/before-base.txt
#   Rscript dev/meme-before.R lane > dev/meme-out/before-lane.txt
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r15",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-meme-lib", libs)
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-meme/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\nfrmtmb.sample",
    format(packageVersion("frmtmb.sample")), "\n")
set.seed(11)
n <- 40
d <- data.frame(x = rnorm(n), z = rnorm(n), sx = 0.3,
                sz = runif(n, 0.2, 0.4), g = factor(rep(1:8, each = 5)),
                h = factor(rep(1:5, 8)))
d$xg <- rnorm(8)[d$g]
d$zg <- rnorm(8)[d$g]
d$sg <- 0.25
d$y <- 1 + 0.5 * d$x - 0.3 * d$z + rnorm(n)
try_show <- function(lab, expr) {
  cat("\n-----", lab, "\n")
  r <- tryCatch(expr, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(r)) print(r)
  invisible(r)
}
cases <- list(
  A = bf(y ~ me(x, sx)),
  B = bf(y ~ me(x, sx) + me(z, sz)),
  E = bf(y ~ me(xg, sg, gr = g) + me(zg, sg, gr = g)),
  F = bf(y ~ me(x, sx) + (1 | h)))
for (nm in names(cases)) {
  f <- cases[[nm]]
  cat("\n==================== case", nm, "\n")
  fit <- frm(f, data = d)
  cat("par_template components:", names(fit$frame$par_template), "\n")
  print(fit$frame$par_template[c("meanme", "logsdme", "thetame")])
  s <- summary(fit)
  try_show("summary()$me", s$me)
  try_show("prior_summary(fit)", as.data.frame(prior_summary(fit)))
  try_show("default_prior(route = 'fit')",
           as.data.frame(default_prior(f, data = d, route = "fit")))
  try_show("default_prior(route = 'sample')",
           as.data.frame(default_prior(f, data = d, route = "sample")))
  try_show("variables(fit) matching me", grep("me", variables(fit),
                                              value = TRUE))
}
try_show("set_prior(class = 'meanme')",
         set_prior("normal(0, 1)", class = "meanme"))
try_show("set_prior(class = 'sdme')",
         set_prior("exponential(1)", class = "sdme"))
try_show("set_prior(class = 'corme')", set_prior("lkj(2)", class = "corme"))
try_show("frm(prior = meanme)",
         frm(cases$A, data = d, prior = set_prior("normal(0, 1)",
                                                  class = "meanme")))

cat("\n==================== frm_sample on case B\n")
fitB <- frm(cases$B, data = d)
dr <- frm_sample(fitB, chains = 1, iter = 300, warmup = 150, seed = 1,
                 refresh = 0)
v <- variables(dr)
cat("variables (non-Xme):", grep("^Xme", v, value = TRUE, invert = TRUE),
    "\n")
cat("Xme names (first 3):", head(grep("^Xme", v, value = TRUE), 3), "\n")
try_show("prior_summary(draws)", as.data.frame(prior_summary(dr)))
try_show("sample_prior = 'yes'", {
  dy <- frm_sample(fitB, chains = 1, iter = 300, warmup = 150, seed = 1,
                   refresh = 0, sample_prior = "yes")
  names(prior_draws(dy))
})
try_show("sample_prior = 'only'", {
  frm_sample(fitB, chains = 1, iter = 300, warmup = 150, seed = 1,
             refresh = 0, sample_prior = "only",
             prior = set_prior("normal(0, 1)", class = "b"))
})
