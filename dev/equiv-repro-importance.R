# Repro for importance = on a plain binomial random-intercept model
# (pair p124, seed 124) and a poisson correlated slope (p069, seed 69):
# frm(importance = ) stops with "did not converge". Run at the survey's
# 200 draws and at 512 and 1000 (the vignettes use 200 to 512), and with
# 20 rounds, to show the stop is not a draw-count artifact.
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("dev/equiv-specs.R", envir = E)
try_fit <- function(lab, expr) {
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  warning:", substr(conditionMessage(w), 1, 150), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(r, "error")) {
    cat(lab, ": ERROR ", substr(conditionMessage(r), 1, 200), "\n", sep = "")
  } else {
    cat(lab, ": logLik ", format(as.numeric(logLik(r)), digits = 10),
        "\n", sep = "")
  }
}
db <- E$eq_dat("binomial", seed = 124)
fb <- function(...) frm(bf(y | trials(ntr) ~ x + (1 | g)), family = binomial(),
                        data = db, ...)
try_fit("binomial Laplace", fb())
try_fit("binomial quadrature", fb(quadrature = TRUE))
for (k in c(200L, 512L, 1000L)) try_fit(paste("binomial importance", k),
                                        fb(importance = k))
try_fit("binomial importance 512, 20 rounds",
        fb(importance = 512L,
           control = frmtmb_control(importance_rounds = 20L)))
dp <- E$eq_dat("poisson", seed = 69)
fp <- function(...) frm(bf(y ~ x + (1 + x | g)), family = poisson(),
                        data = dp, ...)
try_fit("poisson Laplace", fp())
for (k in c(200L, 512L)) try_fit(paste("poisson importance", k),
                                 fp(importance = k))
