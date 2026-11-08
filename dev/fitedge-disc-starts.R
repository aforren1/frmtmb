# Item 2: which starts of the mu body reach -371.54 with disc moved off its
# saddle? Data as dev/fitedge-disc-probe.R (seed 11001). Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
set.seed(11001); n <- 300
dd <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
invisible(factor(rep(1:15, each = 20))); invisible(rnorm(15, 0, 0.5))
lat <- 1.5 * exp(0.6 * dd$x) + 0.7 * dd$z
dd$yo <- ordered(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
fm4 <- bf(yo ~ 1) + nlf(mu ~ b1 * exp(b2 * x)) + nlf(disc ~ exp(d1 * d2 * z)) +
  lf(b1 ~ 1) + lf(b2 ~ 1) + lf(d1 ~ 1) + lf(d2 ~ 1)
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
w <- max(1, diff(range(dd$x)) / 2)
for (b1 in c(0.16, -0.16, 1, -1)) for (b2 in c(0.18 / w, -0.18 / w, 0.5, -0.5)) {
  f <- tryCatch(q(fm4, family = cumulative(), data = dd,
                  start = list(beta = c("b1_(Intercept)" = b1, "b2_(Intercept)" = b2),
                               betad = c(0.3, 0.3))), error = function(e) NULL)
  if (is.null(f)) { cat("b1", b1, "b2", b2, "ERROR\n"); next }
  cat(sprintf("S2 b1 %6.3f b2 %6.3f -> ll %.4f code %d beta %s\n", b1, b2,
              as.numeric(logLik(f)), f$opt$convergence,
              paste(signif(f$estimates$beta, 4), collapse = " ")))
}
# the first fit, before any restart: what it held
ns <- asNamespace("frmtmb")
f0 <- tryCatch(withCallingHandlers(
  frm(fm4, family = cumulative(), data = dd, start = list(betad = c(0, 0)),
      control = frmtmb_control()),
  warning = function(w) { cat("W0:", substr(conditionMessage(w), 1, 200), "\n"); invokeRestart("muffleWarning") },
  message = function(m) { cat("M0:", substr(conditionMessage(m), 1, 200), "\n"); invokeRestart("muffleMessage") }),
  error = function(e) e)
cat("DONE\n")
