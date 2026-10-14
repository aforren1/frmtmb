# Lane postfix, fix round 1: the end of the error each refused scan
# restart gives, on the review's models (dev/postfix-rev-nlscan.R data).
#   Rscript dev/postfix-fix1-msgs.R <lib>
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
tail_of <- function(lbl, expr) {
  m <- tryCatch({ expr; "NO ERROR" }, error = conditionMessage)
  k <- regexpr("A scan of the", m, fixed = TRUE)
  cat(lbl, ":", if (k > 0) substring(m, k) else "(no scan note)", "\n")
}
set.seed(101)
x <- runif(150, 0, 5)
d1 <- data.frame(x = x, y = rpois(150, 20 * exp(-0.5 * x)))
f1 <- bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE)
tail_of("M1c", frm(f1, data = d1, family = poisson("identity"),
                   prior = prior(normal(0, 100), nlpar = "a", lb = 0)))
set.seed(102)
x2 <- rep(c(0.02, 0.06, 0.11, 0.22, 0.56, 1.1), each = 8)
mu2 <- 200 * x2 / (0.1 + x2)
d2 <- data.frame(x = x2, y = rgamma(length(x2), shape = 50, rate = 50 / mu2))
tail_of("M2", frm(bf(y ~ Vm * x / (K + x), Vm + K ~ 1, nl = TRUE),
                  data = d2, family = Gamma("identity")))
set.seed(104)
x4 <- runif(120, 0, 1000)
d4 <- data.frame(x = x4, y = rgamma(120, shape = 20,
                                    rate = 20 / (2e5 * exp(-x4 / 300))))
tail_of("M4", frm(bf(y ~ a * exp(-x / tau), a + tau ~ 1, nl = TRUE),
                  data = d4, family = Gamma("identity")))
set.seed(107)
x7 <- seq(0, 10, length.out = 120)
d7 <- data.frame(x = x7, y = 3 * sin(2 * pi * x7 / 2.5) + rnorm(120, 0, 0.5))
tail_of("M7", frm(bf(y ~ A * sin(6.283185307179586 * x / P), A + P ~ 1,
                     nl = TRUE), data = d7))
set.seed(109)
x9 <- runif(300, 0, 3)
d9 <- data.frame(x = x9)
d9$y <- rgamma(300, shape = 10, rate = 10 / (1 + 2 + 3 * x9 + 0.5 + 1 + 0.2))
tail_of("M9", frm(bf(y ~ p1 + p2 + p3 * x + p4 + p5 + p6,
                     p1 + p2 + p3 + p4 + p5 + p6 ~ 1, nl = TRUE),
                  data = d9, family = Gamma("identity")))
