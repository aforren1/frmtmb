# frmtmb.sample on a cs() FACTOR: does the draws surface carry the new
# bcs_<dummy>[k] names, and do the draws' predictions agree with the
# fit's? The extension reads cs() only through core's cs_offsets_add()
# and brms_coef_table(), with a name-agnostic regex, so the expectation
# is that no change is needed there.
#   Rscript dev/csfactor-sample.R > dev/csfactor-log/sample.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", as.character(packageVersion("frmtmb")),
    "frmtmb.sample", as.character(packageVersion("frmtmb.sample")), "\n")
set.seed(405)
n <- 250
fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(fc)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n)
d <- data.frame(fc = fc,
                yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
go <- function(lab, expr) {
  cat("\n--", lab, "--\n")
  print(tryCatch(expr, error = function(e) paste("ERROR:",
                                                 conditionMessage(e))))
}
ds <- try(frm_sample(bf(yo ~ cs(fc)), family = sratio(), data = d,
                     chains = 2L, iter = 600L, refresh = 0L, seed = 7L))
if (inherits(ds, "try-error")) {
  cat("frm_sample failed\n")
} else {
  go("variables() bcs rows", grep("^bcs", variables(ds), value = TRUE))
  go("fixef(ds) rows", rownames(fixef(ds)))
  go("posterior_epred at the three levels, posterior mean",
     round(apply(posterior_epred(ds, newdata = data.frame(
       fc = factor(c("a", "b", "c")))), c(2L, 3L), mean), 4))
  go("ML fit at the same rows", {
    ff <- frm(bf(yo ~ cs(fc)), family = sratio(), data = d)
    round(frm_linpred(ff, newdata = data.frame(fc = factor(c("a", "b", "c"))),
                      type = "response"), 4)
  })
  go("posterior_epred at the single row factor('c')",
     round(apply(posterior_epred(ds, newdata = data.frame(
       fc = factor("c"))), c(2L, 3L), mean), 4))
}
cat("\ndone\n")
