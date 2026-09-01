# Probe B: the fallback route. A custom_family() whose lpdf calls
# RTMBode::ode() directly, with the per-row times / subject ids / doses
# carried as vreal()/vint() addition terms. This sidesteps the nl
# machinery entirely: every dpar is an ordinary linear predictor.
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib",
            .libPaths()))
library(frmtmb)
library(RTMBode)
source("dev/ode/pk-common.R")

d <- sim_pk()
d$idn <- as.integer(d$id)

pk_family <- custom_family(
  "pk1cmt_oral",
  dpars = c("lka", "lke", "lV", "lsigma"),
  links = list(lka = "identity", lke = "identity", lV = "identity",
               lsigma = "identity"),
  primary_dpars = "lka",
  lpdf = function(y, dpars, aterms) {
    mu <- pk_ode(exp(dpars[["lka"]]), exp(dpars[["lke"]]),
                 exp(dpars[["lV"]]), aterms[["vreal1"]], aterms[["vint1"]],
                 aterms[["vreal2"]])
    RTMB::dnorm(y, mu, exp(dpars[["lsigma"]]), log = TRUE)
  },
  post = list(
    mean_fn = function(dpars, aterms)
      pk_ode(exp(dpars[["lka"]]), exp(dpars[["lke"]]), exp(dpars[["lV"]]),
             aterms[["vreal1"]], aterms[["vint1"]], aterms[["vreal2"]]),
    var_fn = function(dpars, aterms) exp(dpars[["lsigma"]])^2
  ),
  sim = function(dpars, aterms, n)
    rnorm(n, pk_ode(exp(dpars[["lka"]]), exp(dpars[["lke"]]),
                    exp(dpars[["lV"]]), aterms[["vreal1"]],
                    aterms[["vint1"]], aterms[["vreal2"]]),
          exp(dpars[["lsigma"]])),
  type = "continuous")

form <- bf(conc | vreal(time, dose) + vint(idn) ~ 1 + (1 | id),
           lke ~ 1 + (1 | id), lV ~ 1, lsigma ~ 1)

cat("--- fit ---\n")
tt <- system.time(
  fit <- try(frm(form + pk_family, data = d, se = TRUE, verbose = TRUE,
                 start = list(beta = 0, betad = c(log(0.25), log(8), 0))),
             silent = TRUE))
if (inherits(fit, "try-error")) {
  cat("FIT FAILED:\n"); cat(as.character(fit), "\n"); quit(status = 0)
}
cat("elapsed:", tt[["elapsed"]], "s\n")
print(summary(fit))

cat("\nlogLik:", as.numeric(logLik(fit)), "\n")
cat("(nl-route logLik was -60.462931)\n")

cat("\n--- downstream on the custom-family route ---\n")
ck <- function(l, e) {
  r <- try(force(e), silent = TRUE)
  cat(if (inherits(r, "try-error")) "[FAIL] " else "[ok]   ", l,
      if (inherits(r, "try-error")) paste0(": ", sub("\n$", "",
                                                     as.character(r))) else "",
      "\n", sep = "")
  invisible(r)
}
p <- ck("predict()", predict(fit))
if (!inherits(p, "try-error"))
  cat("       max|predict - analytic truth| =",
      format(max(abs(p - d$mu_true)), digits = 3), "\n")
ck("fitted()", fitted(fit))
ck("residuals()", residuals(fit))
ck("simulate()", simulate(fit, nsim = 2))
nd <- d[1:16, ]
ck("predict(newdata=)", predict(fit, newdata = nd))
ck("confint(wald)", confint(fit, method = "wald"))
cat("PROBEB DONE\n")
