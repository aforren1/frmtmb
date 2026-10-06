.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
mode <- "mu"
gen <- function(seed, n = 500) {
  set.seed(seed)
  x <- rnorm(n); cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 2 * x, -1.5 * x) + rlogis(n)
  data.frame(x = x, cls = cls, y = 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5))
}
fam <- mixture(cumulative(), cumulative(), order = mode)
d <- gen(11)
fit <- frm(bf(y ~ x), family = fam, data = d, verbose = FALSE)
print(fit$opt[c("convergence", "message", "objective", "iterations")])
print(fit$opt$par); print(max(abs(fit$obj$gr(fit$opt$par))))
set.seed(1000 + 11)
for (j in 1:11) {
  st <- lapply(fit$frame$par_template, function(v) if (!length(v)) v else v + rnorm(length(v), 0, 1))
  st <- st[c("beta", "betad", grep("^tau_raw", names(st), value = TRUE))]
  f2 <- tryCatch(suppressWarnings(frm(bf(y ~ x), family = fam, data = d, start = st)), error = function(e) NULL)
  if (!is.null(f2)) cat(j, as.numeric(logLik(f2)), round(f2$opt$par, 3), f2$opt$convergence, "\n")
}
