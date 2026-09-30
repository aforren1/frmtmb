# Is a flexible ordinal fit without a modeled disc the same fit, to the
# bit, on the lane and on the base build? Each arm writes its estimates,
# logLik and fitted probabilities; arm "compare" reads both. Seed
# 20260930. Usage: Rscript dev/ordinal-bitwise.R lane|base|compare
arm <- commandArgs(TRUE)[1]
out <- function(a) sprintf("dev/ordinal-bitwise-%s.rds", a)
if (identical(arm, "compare")) {
  a <- readRDS(out("base"))
  b <- readRDS(out("lane"))
  for (nm in names(a)) {
    cat(sprintf("%-26s par %s  logLik %s  fitted %s  osa %s\n", nm,
                identical(a[[nm]]$par, b[[nm]]$par),
                identical(a[[nm]]$ll, b[[nm]]$ll),
                identical(a[[nm]]$P, b[[nm]]$P),
                identical(a[[nm]]$osa, b[[nm]]$osa)))
  }
  quit(save = "no")
}
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib", libs)
}
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
set.seed(20260930)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
fits <- list(
  cumulative = frm(y ~ x, family = cumulative(), data = d),
  cumulative_probit = frm(y ~ x, family = cumulative("probit"), data = d),
  sratio_cs = frm(y ~ x + cs(z), family = sratio(), data = d),
  cratio_cloglog = frm(y ~ x, family = cratio("cloglog"), data = d),
  acat_cs = frm(y ~ x + cs(z), family = acat(), data = d),
  cumulative_thres_gr = frm(y | thres(gr = g) ~ x, family = cumulative(),
                            data = d),
  hurdle = frm(yh ~ x, family = hurdle_cumulative(), data = d)
)
res <- lapply(fits, function(f) {
  list(par = f$opt$par, ll = as.numeric(logLik(f)),
       P = fitted(f)[, "Estimate", ],
       osa = tryCatch({
         set.seed(1)
         residuals(f, type = "osa")[, "Estimate"]
       }, error = function(e) conditionMessage(e)))
})
saveRDS(res, out(arm))
cat("wrote", out(arm), "\n")
