# Lane prilabel: which standard deviations frm_sample()'s defaults keep
# when the user writes a narrower class "sd" prior. brms 2.23.0 keeps
# its student_t(3, 0, 2.5) default on every sd the user did not name
# (dev/prilabel-log/brms-ref.txt: student_t on sd_1[1] beside a
# per-coefficient normal on sd_1[2]). Arg 1: "lane" or "base".
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", dirname(find.package("frmtmb")), "; frmtmb.sample",
    dirname(find.package("frmtmb.sample")), "\n")
set.seed(3)
n <- 200
dd <- data.frame(x = rnorm(n), g = factor(sample(1:12, n, TRUE)),
                 h = factor(sample(1:6, n, TRUE)))
dd$y <- 1 + dd$x + rnorm(12)[dd$g] + rnorm(6)[dd$h] + rnorm(n)
fit <- frm(bf(y ~ x + (1 + x | g) + (1 | h)), family = gaussian(), data = dd)
show <- function(pr, lab) {
  rp <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, pr))
  e <- Filter(function(e) identical(e$comp, "theta") &&
                !identical(e$dist$kind, "lkj"), rp$ri$entries)
  d <- vapply(e, function(e) paste0(e$dist$kind, "(", e$dist$scale, ")"),
              "")
  i <- vapply(e, function(e) as.integer(e$idx), 0L)
  o <- order(i)
  cat(sprintf("%-34s theta %s\n", lab,
              paste(paste0(i[o], "=", d[o]), collapse = "  ")))
}
cat("theta 1, 2 = sd_g__Intercept, sd_g__x; 4 = sd_h__Intercept\n")
show(NULL, "no user prior")
r <- tryCatch(show(set_prior("normal(0, 0.5)", class = "sd", group = "g",
                             coef = "x"), "sd, group g, coef x"),
              error = function(e) cat("sd, group g, coef x: ERROR",
                                      conditionMessage(e), "\n"))
show(set_prior("normal(0, 0.5)", class = "sd", group = "g"), "sd, group g")
show(set_prior("normal(0, 0.5)", class = "sd"), "sd (class-wide)")
