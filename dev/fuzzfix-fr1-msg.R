# Fix round 1 copy of dev/fuzzfix-rev-msg.R on the lane library: the full text of the two new refusals and of the crossed
# quadrature refusal, as a user sees them.
.libPaths(c("C:/Users/adf44/source/r/wt-fuzzfix-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(20261021)
n <- 72
d <- data.frame(g = factor(rep(1:12, each = 6)), tt = factor(rep(1:6, 12)),
                x = rnorm(n), w = runif(n, 0.5, 2))
d$obs <- factor(seq_len(n))
d$y <- rgamma(n, shape = 3, rate = 3 / exp(1 + 0.4 * d$x))
d$ly <- log(d$y)
p <- function(e) cat(conditionMessage(e), "\n\n")
tryCatch(frm(bf(y | weights(w) ~ x + ar1(tt + 0 | g)),
             family = Gamma("log"), data = d), error = p)
tryCatch(suppressWarnings(frm(bf(ly | weights(w) ~ x + (1 | obs)),
                              family = gaussian(), data = d)), error = p)
tryCatch(frm(bf(y | weights(w) ~ x + (1 | obs)), family = Beta(),
             data = transform(d, y = plogis(ly))), error = p)
