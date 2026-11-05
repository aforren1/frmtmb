# The sweep's beta_binomial row read logLik 6.6e12 on the ridge-restart rule: what
# did the restart reach? The sweep's own draws (dev/equiv-nl-sweep.R, seed
# 20261029), replayed up to that family.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- readLines("dev/equiv-nl-sweep.R")
a <- grep("^ft <- frm_compat_features", src)
b <- grep("^lhs <- function", src)
eval(parse(text = src[a:(b - 1L)]))
for (f in fams) {
  y <- ydat(f)
  if (identical(f, "beta_binomial")) break
}
d <- data.frame(x = x, ntr = 10L)
d$y <- y
body <- bf(y | trials(ntr) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)
w <- character(0); m <- character(0)
fit <- withCallingHandlers(frm(body, family = beta_binomial(), data = d),
  warning = function(c) { w <<- c(w, substr(conditionMessage(c), 1, 200)); invokeRestart("muffleWarning") },
  message = function(c) { m <<- c(m, substr(conditionMessage(c), 1, 200)); invokeRestart("muffleMessage") })
cat("logLik", as.numeric(logLik(fit)), "conv", fit$opt$convergence, "\n")
print(unlist(fit$estimates[c("beta", "betad")]))
cat("w:", w, "\nm:", m, "\n")
lin <- frm(bf(y | trials(ntr) ~ x), family = beta_binomial(), data = d)
cat("linear", as.numeric(logLik(lin)), "conv", lin$opt$convergence, "\n")
print(unlist(lin$estimates[c("beta", "betad")]))
