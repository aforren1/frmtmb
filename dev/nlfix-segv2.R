# A saddle-table case that dies with exit 139: which step?
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(9202)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
set.seed(1)
d$yh <- ordered(ifelse(runif(n) < 0.2, 0L, as.integer(d$y)), levels = 0:4)
set.seed(9205)
d$yg <- 1.5 * exp(0.6 * d$x) + rnorm(n, 0, 0.5)
d$yc <- rpois(n, 1.5 * exp(0.6 * d$x))
f <- bf(yc ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)
cat("default\n"); flush.console()
r <- tryCatch(suppressWarnings(suppressMessages(frm(f, family = poisson("identity"), data = d))), error = function(e) e)
cat(if (inherits(r, "error")) conditionMessage(r) else as.numeric(logLik(r)), "\n"); flush.console()
cat("started\n"); flush.console()
r <- tryCatch(suppressWarnings(frm(f, family = poisson("identity"), data = d, start = list(beta = c("b1_(Intercept)" = 0.3, "b2_(Intercept)" = 0.3)))), error = function(e) e)
cat(if (inherits(r, "error")) conditionMessage(r) else as.numeric(logLik(r)), "\n"); flush.console()
for (s in list(c(.16, -.18), c(-.16, .18), c(-.16, -.18), c(.16, .18))) {
  cat("start", s, "\n"); flush.console()
  r <- tryCatch(suppressWarnings(frm(f, family = poisson("identity"), data = d, start = list(beta = c("b1_(Intercept)" = s[1], "b2_(Intercept)" = s[2])))), error = function(e) e)
  cat(if (inherits(r, "error")) substr(conditionMessage(r), 1, 100) else as.numeric(logLik(r)), "\n"); flush.console()
}
