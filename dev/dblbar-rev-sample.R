# Reviewer, lane dblbar: frmtmb.sample on (1 + f || g) against the same
# model as three blocks over numeric dummies. Names, priors, draws,
# hypothesis(), posterior_predict() at new levels, as_tmbstan().
# Data seed 20261013 (dev/dblbar-data.R); sampler seed 11, 4 chains of
# 1500. Output: dev/dblbar-rev-log/sample.txt
.libPaths(c("C:/Users/adf44/source/r/wt-dblbar-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
sink(file.path(wt, "dev/dblbar-rev-log/sample.txt"), split = TRUE)
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
dd <- dblbar_data()
dd$fb <- as.numeric(dd$f == "b")
dd$fc <- as.numeric(dd$f == "c")
tc <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
A <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
B <- frm(bf(y ~ f + (1 | g) + (0 + fb | g) + (0 + fc | g)),
         family = gaussian(), data = dd)
pA <- as.data.frame(default_prior(bf(y ~ f + (1 + f || g)), data = dd,
                                  family = gaussian(), route = "sample"))
print(pA[pA$class %in% c("sd", "cor"), c("prior", "class", "coef", "group")])
cat("coef sd prior:", substr(tc(frm(bf(y ~ f + (1 + f || g)),
                                    family = gaussian(), data = dd,
                                    prior = set_prior("normal(0, 1)",
                                                      class = "sd",
                                                      coef = "fb",
                                                      group = "g"))), 1, 90),
    "\n")
dA <- frm_sample(A, chains = 4, iter = 1500, seed = 11, refresh = 0)
dB <- frm_sample(B, chains = 4, iter = 1500, seed = 11, refresh = 0)
vA <- VarCorr(dA)
cat("VarCorr(dA)$g elements:", names(vA$g), "rows",
    rownames(vA$g$sd), "\n")
print(vA$g$sd)
vB <- VarCorr(dB)
print(vB$g$sd)
cat("draws names with cor:", grep("cor", variables(dA), value = TRUE), "\n")
cat("hypothesis sd_g__fb > sd_g__fc:\n")
h <- tc(hypothesis(dA, "sd_g__fb > sd_g__fc", class = NULL))
print(if (is.character(h)) h else h$hypothesis)
cat("hypothesis on a cor_:", substr(tc(hypothesis(dA,
  "cor_g__Intercept__fb > 0", class = NULL))[1], 1, 120), "\n")
nd <- dd[c(1, 5, 9, 40), ]
nd$g <- c("NEW1", "NEW2", "NEW1", "NEW3")
for (snl in c("gaussian", "old_levels", "uncertainty")) {
  set.seed(7)
  a <- tc(posterior_epred(dA, newdata = nd, allow_new_levels = TRUE,
                          sample_new_levels = snl))
  set.seed(7)
  b <- tc(posterior_epred(dB, newdata = nd, allow_new_levels = TRUE,
                          sample_new_levels = snl))
  if (is.character(a) || is.character(b)) {
    cat(snl, "A:", substr(a[1], 1, 100), "| B:", substr(b[1], 1, 100), "\n")
    next
  }
  cat(sprintf("epred %-11s A mean %s sd %s\n", snl,
              paste(sprintf("%.3f", colMeans(a)), collapse = " "),
              paste(sprintf("%.3f", apply(a, 2, sd)), collapse = " ")))
  cat(sprintf("epred %-11s B mean %s sd %s\n", snl,
              paste(sprintf("%.3f", colMeans(b)), collapse = " "),
              paste(sprintf("%.3f", apply(b, 2, sd)), collapse = " ")))
}
set.seed(8)
pp <- tc(posterior_predict(dA, newdata = nd, allow_new_levels = TRUE))
cat("posterior_predict dim:", if (is.character(pp)) pp else dim(pp),
    " means", if (!is.character(pp)) sprintf("%.3f", colMeans(pp)), "\n")
ts <- tc(as_tmbstan(A, chains = 1, iter = 400, seed = 1, refresh = 0))
cat("as_tmbstan:", if (is.character(ts)) ts else class(ts), "\n")
if (!is.character(ts)) {
  cat("as_tmbstan pars:", head(names(ts), 12), "\n")
}
sink()
