# Reviewer probe: mi() refusal, GEV mu on a log link with negative data,
# discrete_weibull identity link outside (0, 1) on the tape, and how a
# mixture reports a scalar xi.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
try1 <- function(lab, expr) {
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  [warning]", substr(conditionMessage(w), 1, 160), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 300)))
  if (inherits(r, "frmtmb")) r <- sprintf("fit: conv %d logLik %.6f",
                                          r$opt$convergence, as.numeric(logLik(r)))
  cat(lab, ":", if (is.character(r)) r else class(r)[1], "\n")
}
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
set.seed(1301)
n <- 200
d <- data.frame(x = rnorm(n))
d$yf <- exp(1 + 0.3 * d$x) / gamma(1 - 1 / 4) * rexp(n)^(-1 / 4)
d$yg <- rg(n, -2 + 0.5 * d$x, 0.8, 0.1)
d$yf[1:5] <- NA
try1("frechet with mi()", frm(yf | mi() ~ x, family = frechet(), data = d))
try1("GEV link = 'log', mean(y) < 0",
     frm(yg ~ x, family = gen_extreme_value(link = "log"), data = d))
d$yg2 <- d$yg + 6
try1("GEV link = 'log', mean(y) > 0",
     frm(yg2 ~ x, family = gen_extreme_value(link = "log"), data = d))
fam <- discrete_weibull(link = "identity")
f <- function(p) fam$lpdf(c(0, 2), list(mu = p[1], shape = 1), list())
tp <- RTMB::MakeTape(f, 0.5)
cat("dw identity link on the tape, mu = 1.2:", tp(1.2), "; mu = -0.2:", tp(-0.2),
    "\n")
d$yc <- c(rg(n / 2, 0, 1, 0.1), rg(n / 2, 6, 1, 0.1))
fm <- frm(yc ~ 1, family = mixture(gen_extreme_value(), gen_extreme_value()),
          data = d)
print(summary(fm))
