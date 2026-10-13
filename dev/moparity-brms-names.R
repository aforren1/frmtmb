# brms 2.23.0's names for mo() interaction columns, by written order.
# Usage: Rscript dev/moparity-brms-names.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(brms); library(frmtmb)})
set.seed(1)
n <- 120
dd <- data.frame(x = sample(0:3, n, TRUE), z = rnorm(n), w = rnorm(n),
                 f2 = factor(sample(c("a", "b"), n, TRUE)))
dd$y <- rnorm(n)
for (fo in list(y ~ z + mo(x):z, y ~ mo(x):z, y ~ mo(x) * z, y ~ z * mo(x),
                y ~ w + mo(x):z + mo(x), y ~ f2 + mo(x):f2,
                y ~ mo(x):z:w, y ~ mo(x) * z * w)) {
  cat("\n==", deparse1(fo), "\n")
  dp <- brms::default_prior(fo, data = dd)
  cat("brms b  :", dp$coef[dp$class == "b" & nzchar(dp$coef)], "\n")
  cat("brms simo:", dp$coef[dp$class == "simo"], "\n")
  sc <- brms::stancode(fo, data = dd)
  cat("brms eta:", grep("mu\\[n\\] \\+=", strsplit(sc, "\n")[[1]],
                        value = TRUE), "\n")
  if (!grepl("f2", deparse1(fo))) r <- tryCatch({
    f <- frm(bf(fo), family = gaussian(), data = dd)
    paste(names(fixef(f)[, 1]), collapse = " ")
  }, error = function(e) paste("ERROR", conditionMessage(e)))
  cat("frmtmb  :", r, "\n")
}
