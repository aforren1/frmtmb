# The saddle's restarts one by one: at 3e4 rows on the data of dev/nlfix-cost2.R,
# and at 300 rows on a body whose slope has the other sign. Which restarts does a
# stop at two agreeing ones save, and what would they reach?
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")), "\n")
ns <- asNamespace("frmtmb")
ctl <- frmtmb_control(optCtrl = list(iter.max = 0))
each <- function(lab, f, fam, d) {
  cat("---", lab, "\n")
  at0 <- suppressWarnings(frm(f, family = fam, data = d, control = ctl))
  sd <- ns$nl_saddle_find(at0, use_curvature = FALSE)$saddle
  tm <- system.time(fit1 <- suppressWarnings(suppressMessages(frm(f, family = fam, data = d))))[["elapsed"]]
  cat(sprintf("frm() whole: %.2f s ll %.4f code %d\n", tm, as.numeric(logLik(fit1)), fit1$opt$convergence))
  for (st in ns$nl_saddle_starts(at0$frame, sd, NULL)) {
    w <- character(0)
    tm <- system.time(r <- withCallingHandlers(suppressMessages(
      frm(f, family = fam, data = d, start = st, control = frmtmb_control())),
      warning = function(c) { w <<- c(w, substr(conditionMessage(c), 1, 30)); invokeRestart("muffleWarning") }))[["elapsed"]]
    cat(sprintf("start %s: %.2f s ll %.4f code %d beta %s flat %s warn[%s]\n",
                paste(round(unlist(st$beta), 3), collapse = ","), tm, as.numeric(logLik(r)),
                r$opt$convergence, paste(round(r$estimates$beta, 3), collapse = ","),
                any(grepl("Standard errors", w)), paste(unique(w), collapse = " | ")))
  }
}
set.seed(6401)
n <- 3e4
d <- data.frame(x = rnorm(n))
d$yo <- cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
each("linear truth, b1 * exp(b2 * x), 3e4", bf(yo ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), cumulative(), d)
set.seed(11001)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
d$y2 <- ordered(cut(1.2 * exp(-0.5 * d$z) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
each("negative slope, c1 * exp(c2 * z), 300", bf(y2 ~ c1 * exp(c2 * z), c1 + c2 ~ 1, nl = TRUE), cumulative(), d)
each("the same under acat", bf(y2 ~ c1 * exp(c2 * z), c1 + c2 ~ 1, nl = TRUE), acat(), d)
d$y3 <- 1 + 0.5 * d$x + rnorm(n)
each("a + b1 * b2 * x, 300", bf(y3 ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE), gaussian(), d)
