# Lane dblbar: are the nested-group order and the prior-table group
# labels of g1/g2 and mm() a property of `|` terms too? Base build only.
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(frmtmb))
dd <- dblbar_data()
sink(file.path(wt, "dev/dblbar-log/preexist.txt"), split = TRUE)
for (form in list(y ~ x + (1 + x | g1/g2), y ~ x + (1 + x | mm(m1, m2)))) {
  cat("\n==", deparse1(form), "\n")
  fit <- frm(bf(form), family = gaussian(), data = dd)
  print(grep("^sd_|^cor_", variables(fit), value = TRUE))
  p <- as.data.frame(default_prior(bf(form), data = dd, family = gaussian()))
  print(p[p$class %in% c("sd", "cor"), c("class", "coef", "group")])
  r <- tryCatch({
    frm(bf(form), family = gaussian(), data = dd,
        prior = set_prior("normal(0, 1)", class = "sd",
                          group = if (grepl("mm", deparse1(form))) "mmm1m2"
                                  else "g1:g2"))
    "brms's group label accepted"
  }, error = function(e) conditionMessage(e))
  cat("set_prior with brms's group label:", r, "\n")
}
sink()
