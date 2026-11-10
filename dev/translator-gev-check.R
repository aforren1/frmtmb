# Pair p021 (gen_extreme_value x mixture): is frm()'s log-likelihood at its
# estimate the mixture's own, computed here by hand with a GEV density
# that is zero outside its support, and where does each component's xi
# sit against the interval brms's scale_xi() allows?
#
#   Rscript dev/translator-gev-check.R
#
# From the worktree root; reads FRMTMB_EQUIV_LIB as dev/equiv-measure.R.
lib <- Sys.getenv("FRMTMB_EQUIV_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
eq_env <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = eq_env)
}
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
spec <- eq_env$eq_specs[["gen_extreme_value x mixture"]]()
fit <- eq_env$eq_frm_fit(spec)
str(fit$estimates)
dp <- frmtmb:::eval_dpars(fit)[[1L]]
y <- spec$data$y

# GEV log density with location m, scale s, shape xi; -Inf off the support
lgev <- function(y, m, s, xi) {
  z <- (y - m) / s
  t <- 1 + xi * z
  out <- rep(-Inf, length(y))
  ok <- t > 0
  out[ok] <- -log(s) - (1 + 1 / xi) * log(t[ok]) - t[ok]^(-1 / xi)
  out
}
lp <- function(k) {
  # sigma and xi have no predictor here, so each is one number
  lgev(y, as.numeric(dp[[paste0("mu", k)]]),
       unique(as.numeric(dp[[paste0("sigma", k)]])),
       unique(as.numeric(dp[[paste0("xi", k)]])))
}
th <- frmtmb:::fixef_by_dpar(fit)[["theta1"]][["(Intercept)"]]
w <- exp(c(th, 0)) / sum(exp(c(th, 0)))
l1 <- lp(1) + log(w[1])
l2 <- lp(2) + log(w[2])
m <- pmax(l1, l2)
hand <- sum(m + log(exp(l1 - m) + exp(l2 - m)))
cat(sprintf("logLik(fit) %.10f  by hand %.10f  diff %.3g\n",
            as.numeric(logLik(fit)), hand, as.numeric(logLik(fit)) - hand))
for (k in 1:2) {
  x <- (y - as.numeric(dp[[paste0("mu", k)]])) /
    as.numeric(dp[[paste0("sigma", k)]])
  b <- c(-1 / min(x), -1 / max(x))
  xi <- unique(as.numeric(dp[[paste0("xi", k)]]))
  cat(sprintf(paste0("component %d: xi %.6f, brms interval (%.6f, %.6f); ",
                     "rows off its support %d of %d\n"), k, xi, min(b),
              max(b), sum(1 + xi * x <= 0), length(y)))
}
