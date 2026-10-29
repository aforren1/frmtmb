# Review of lane meme: follow-ups to dev/review-meme-mach.R.
#   Rscript dev/review-meme-mach2.R > dev/review-meme-out/mach2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-meme/dev/review-meme-shapes.R")
d <- rv_data()
q <- function(expr) suppressMessages(suppressWarnings(expr))
try_ <- function(lab, expr) {
  r <- tryCatch(q(expr), error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("--", lab, "\n")
  if (is.character(r) && length(r) == 1) cat("  ", r, "\n") else print(r)
  invisible(r)
}
f_me <- bf(y ~ me(x, sx) + me(z, sz) + w)
cat("#### duplicates\n")
try_("frm: sdme class row twice", {
  f <- frm(f_me, data = d, prior = set_prior("exponential(1)", class = "sdme") +
             set_prior("exponential(5)", class = "sdme"))
  as.data.frame(prior_summary(f))[, c("prior", "class", "coef")]})
try_("frm: corme class row twice", {
  f <- frm(f_me, data = d, prior = set_prior("lkj(2)", class = "corme") +
             set_prior("lkj(5)", class = "corme"))
  as.data.frame(prior_summary(f))[, c("prior", "class", "coef")]})
try_("frm: b class row twice", {
  f <- frm(f_me, data = d, prior = set_prior("normal(0, 1)", class = "b") +
             set_prior("normal(0, 5)", class = "b"))
  as.data.frame(prior_summary(f))[, c("prior", "class", "coef")]})
try_("frm: sigma twice", {
  f <- frm(f_me, data = d, prior = set_prior("exponential(1)", class = "sigma") +
             set_prior("exponential(5)", class = "sigma"))
  as.data.frame(prior_summary(f))[, c("prior", "class", "coef")]})
try_("check_prior_slots(sdme twice)",
     check_prior_slots(set_prior("exponential(1)", class = "sdme") +
                         set_prior("exponential(5)", class = "sdme")))

cat("\n#### prior draws under a coef bound: class b, for comparison\n")
ds <- q(frm_sample(f_me, data = d, family = gaussian(),
                   prior = set_prior("normal(0, 1)", class = "b") +
                     set_prior("normal(0, 1)", class = "b", coef = "w",
                               ub = 0.1),
                   sample_prior = "yes", chains = 2, iter = 1000, seed = 3,
                   refresh = 0))
pd <- prior_draws(ds)
cat("  prior draw names:", names(pd), "\n")
for (nm in grep("^b", names(pd), value = TRUE)) {
  cat(sprintf("  prior %s: max %.3f  n = %d\n", nm, max(pd[[nm]]),
              length(pd[[nm]])))
}
cat(sprintf("  posterior max b_w %.3f\n", max(as.matrix(ds)[, "b_w"])))
cat("\n#### prior draws under a class bound on sdme and meanme\n")
ds3 <- q(frm_sample(f_me, data = d, family = gaussian(),
                    prior = set_prior("exponential(1)", class = "sdme",
                                      ub = 0.5) +
                      set_prior("normal(0, 1)", class = "meanme", lb = 0),
                    sample_prior = "yes", chains = 2, iter = 1000, seed = 3,
                    refresh = 0))
pd3 <- prior_draws(ds3)
cat("  names:", names(pd3), "\n")
cat(sprintf("  prior sdme_1 max %.3f; prior meanme_1 min %.3f\n",
            max(pd3$sdme_1), min(pd3$meanme_1)))
m3 <- as.matrix(ds3)
cat(sprintf("  posterior sdme max %.3f; meanme min %.3f\n",
            max(m3[, c("sdme_mex", "sdme_mez")]),
            min(m3[, c("meanme_mex", "meanme_mez")])))

cat("\n#### hypothesis(class = NULL) on meanme\n")
ds2 <- q(frm_sample(f_me, data = d, family = gaussian(),
                    prior = set_prior("normal(0, 1)", class = "meanme") +
                      set_prior("exponential(1)", class = "sdme"),
                    sample_prior = "yes", chains = 2, iter = 1000, seed = 3,
                    refresh = 0))
for (h in c("meanme_mex = 0", "meanme_mex > 0", "sdme_mex = 1",
            "corme__mex__mez = 0", "sigma = 1")) {
  try_(h, hypothesis(ds2, h, class = NULL)$hypothesis[, c("Estimate",
                                                         "Evid.Ratio",
                                                         "Post.Prob")])
}
cat("DONE\n")
