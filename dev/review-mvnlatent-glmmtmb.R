# us_sigma()'s correlation order against glmmTMB's get_cor(), d = 3, 4.
.libPaths(c("C:/Users/adf44/source/r/wt-mvnlatent-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("glmmTMB", format(packageVersion("glmmTMB")), "\n")
print(body(glmmTMB::get_cor))
for (d in 3:4) {
  set.seed(d)
  tc <- round(rnorm(d * (d - 1) / 2), 2)
  S <- as.matrix(us_sigma(c(numeric(d), tc), d))
  G <- glmmTMB::get_cor(tc, return_val = "mat")
  cat("d =", d, " equal to glmmTMB get_cor:", isTRUE(all.equal(S, G)), "\n")
}
