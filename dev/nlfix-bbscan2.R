# The same data on the base build: is a positive logLik of beta_binomial() with
# no overdispersion older than this lane, and on the linear spelling too?
lib <- Sys.getenv("NLFIX_LIB", "C:/Users/adf44/source/r/rellib-r16")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat(find.package("frmtmb"), "\n")
for (s in c(2, 4, 6, 9, 12)) {
  set.seed(s)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), ntr = 10L)
  d$y <- rbinom(n, 10, plogis(0.3 + 0.5 * d$x))
  a <- suppressWarnings(frm(bf(y | trials(ntr) ~ x), family = beta_binomial(), data = d))
  b <- suppressWarnings(frm(bf(y | trials(ntr) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                            family = beta_binomial(), data = d))
  cat(s, "linear", as.numeric(logLik(a)), a$estimates$betad, "nl", as.numeric(logLik(b)), b$estimates$betad, "\n")
}
