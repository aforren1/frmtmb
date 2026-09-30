# Reviewer, claim 4: rate() with cens() on poisson. brms_lp_check found
# a gap of 32.66 (dev/aterms2-rev-log-03-brms.txt). Here: the IDENTITY
# rate(time) + cens(cc) against offset(log(time)) + cens(cc) under the
# log link, at the same parameter point, and brms's Stan code.
# Seed 306 data. Log: dev/aterms2-rev-log-03c-ratecens.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(306)
n <- 250
dr <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4), wt = runif(n, 0.5, 2))
dr$y <- rpois(n, exp(0.4 + 0.3 * dr$x) * dr$time)
dr$cc <- rep(c(0, 0, 1, -1), length.out = n)
fr <- q(frm(y | rate(time) + cens(cc) ~ x, data = dr, family = poisson()))
fo <- q(frm(y | cens(cc) ~ x + offset(log(time)), data = dr, family = poisson()))
p <- fo$opt$par
cat("fn rate+cens", sprintf("%.10f", fr$obj$fn(p)), " fn offset+cens",
    sprintf("%.10f", fo$obj$fn(p)), "\n")
cat("logLik rate+cens", format(logLik(fr), digits = 12), " offset+cens",
    format(logLik(fo), digits = 12), "\n")
print(rbind(rate = fixef(fr)[, 1], offset = fixef(fo)[, 1]))
# by hand at the offset fit's estimates, frmtmb's cens convention read
# from the uncensored model: left (-1) is P(Y <= y), right (1) P(Y >= y)
b <- fixef(fo)[, 1]
mu <- exp(b[1] + b[2] * dr$x) * dr$time
ll <- ifelse(dr$cc == 0, dpois(dr$y, mu, log = TRUE),
       ifelse(dr$cc == -1, ppois(dr$y, mu, log.p = TRUE),
              ppois(dr$y - 1, mu, lower.tail = FALSE, log.p = TRUE)))
cat("by hand, P(Y >= y) right:", sum(ll), "\n")
ll2 <- ifelse(dr$cc == 0, dpois(dr$y, mu, log = TRUE),
        ifelse(dr$cc == -1, ppois(dr$y, mu, log.p = TRUE),
               ppois(dr$y, mu, lower.tail = FALSE, log.p = TRUE)))
cat("by hand, P(Y > y) right:", sum(ll2), "\n")
# the same with mu WITHOUT the exposure on the censored rows
mu0 <- exp(b[1] + b[2] * dr$x)
ll3 <- ifelse(dr$cc == 0, dpois(dr$y, mu, log = TRUE),
        ifelse(dr$cc == -1, ppois(dr$y, mu0, log.p = TRUE),
               ppois(dr$y - 1, mu0, lower.tail = FALSE, log.p = TRUE)))
cat("by hand, exposure dropped on censored rows:", sum(ll3), "\n")
cat("-fn(rate+cens) at the same point:", -fr$obj$fn(p), "\n")
code <- brms::stancode(brms::bf(y | rate(time) + cens(cc) ~ x), data = dr,
                       family = poisson())
cl <- strsplit(code, "\n")[[1]]
cat("\nbrms Stan lines with target:\n")
cat(cl[grepl("target|denom|lccdf|lcdf", cl)], sep = "\n")
code0 <- brms::stancode(brms::bf(y | cens(cc) ~ x + offset(log(time))),
                        data = dr, family = poisson())
cl0 <- strsplit(code0, "\n")[[1]]
cat("\nbrms offset + cens lines with target:\n")
cat(cl0[grepl("target|lccdf|lcdf", cl0)], sep = "\n")
