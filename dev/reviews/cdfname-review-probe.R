# Review probe: the `$` reader, frmtmb_family()'s `...`, and routes the
# lane's test does not list.
#   Rscript dev/reviews/cdfname-review-probe.R
.libPaths(c("C:/Users/adf44/source/r/wt-cdfname-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
say <- function(label, expr) {
  r <- tryCatch({
    v <- suppressMessages(force(expr))
    if (is.null(v)) "value NULL" else paste("value of class", class(v)[1L])
  }, error = function(e) {
    sprintf("ERROR [%s] %s", class(e)[1L],
            substr(gsub("\n", " ", conditionMessage(e)), 1, 150))
  })
  cat(sprintf("P %-44s %s\n", label, r))
}
pois_family <- function(...) {
  frmtmb_family(
    "pois_rev", dpars = "mu", links = list(mu = "log"),
    type = "discrete", accepts_aterms = c("weights", "cens", "trunc"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dpois(y, dpars[["mu"]], log = TRUE)
    },
    init_dpars = list(mu = function(y, aterms) mean(y)), ...)
}
pcdf <- function(q, dpars, aterms) RTMB::ppois(q, dpars[["mu"]])
old <- function() {
  f <- pois_family()
  f[["lcdf"]] <- pcdf
  f
}
set.seed(3088)
d <- data.frame(y = rpois(80, 4), y2 = rpois(80, 3))
d$cc <- ifelse(d$y > 6, "right", "none")
d$y <- pmin(d$y, 6)

cat("-- item 2: `$` and the slot names\n")
f0 <- pois_family()
cat("P names hold cdf with NULL value:", "cdf" %in% names(f0),
    "| cdf_log:", "cdf_log" %in% names(f0), "\n")
say("f0$cdf (no CDF given)", f0$cdf)
fd <- as_frmtmb_family(poisson())
fd[["cdf"]] <- NULL
say("poisson() with cdf element dropped, $cdf", fd$cdf)
fl <- as_frmtmb_family(poisson())
fl[["cdf_log"]] <- NULL
say("poisson() with cdf_log dropped, $cdf_log", fl$cdf_log)
say("$cd (prefix of two)", fd$cd)
say("old()$lcdf (stored element)", old()$lcdf)

cat("-- item 4: frmtmb_family()'s `...`\n")
say("cdff = (misspelled)", pois_family(cdff = pcdf))
say("lpdff = (misspelled)", pois_family(lpdff = pcdf))
say("lccd = (prefix of lccdf)", pois_family(lccd = function(q, dpars, aterms) 0))
say("cdf_l = (prefix of cdf_log)", pois_family(cdf_l = pcdf))
say("lcdf = NULL", pois_family(lcdf = NULL))
say("custom_family(lcdf = )", custom_family(
  "x", dpars = "mu", links = list(mu = "log"),
  lpdf = function(y, dpars, aterms) 0, lcdf = pcdf))
say("do.call with lcdf", do.call(frmtmb_family, list(
  "x", dpars = "mu", links = list(mu = "log"),
  lpdf = function(y, dpars, aterms) 0, lcdf = pcdf)))
say("one extra positional", do.call(frmtmb_family, c(list("x"),
  rep(list(NULL), length(formals(frmtmb_family)) - 1L), list(3))))

cat("-- item 3: more routes\n")
say("mvbf(bf(old), bf(poisson))", frm(mvbf(bf(y ~ 1, family = old()),
                                          bf(y2 ~ 1, family = poisson())),
                                     data = d))
mv <- mvbf(bf(y ~ 1, family = pois_family(cdf = pcdf)),
           bf(y2 ~ 1, family = poisson()))
say("mvbf edited after as_frmtmb_family", {
  mv$forms[[1]]$family[["lcdf"]] <- pcdf
  frm(mv, data = d)
})
say("mixture of old inside mvbf", mvbf(bf(y ~ 1) + mixture(old(), poisson()),
                                       bf(y2 ~ 1, family = poisson())))
fit <- frm(y | cens(cc) ~ 1, family = pois_family(cdf = pcdf), data = d)
say("update(fit, family = old())", update(fit, family = old()))
fit_old <- fit
for (p in list(c("bform", "family"), c("spec", "responses", "y", "family"),
               c("frame", "spec", "responses", "y", "family"))) {
  fam <- fit_old[[p]]
  fam[["lcdf"]] <- fam[["cdf"]]
  fam[["cdf"]] <- NULL
  fit_old[[p]] <- fam
}
say("old-shaped fit: residuals(type = osa)", residuals(fit_old, type = "osa"))
say("old-shaped fit: predict(newdata)", predict(fit_old, newdata = d[1:3, ]))
say("old-shaped fit: update()", update(fit_old))
if (requireNamespace("frmtmb.sample", quietly = TRUE)) {
  cat("frmtmb.sample:", find.package("frmtmb.sample"), "\n")
  fs <- getNamespaceExports("frmtmb.sample")
  cat("P frmtmb.sample exports frm_sample:", "frm_sample" %in% fs, "\n")
  if ("frm_sample" %in% fs) {
    say("frm_sample(family = old())", frmtmb.sample::frm_sample(
      y | cens(cc) ~ 1, family = old(), data = d, chains = 1, iter = 20))
  }
}
