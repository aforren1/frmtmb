# Reviewer: the student-t latent weighted case on rellib-r12 (0.74.0).
.libPaths(c("C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(20261020)
n <- 72
d <- data.frame(x = rnorm(n)); d$obs <- factor(seq_len(n))
eta <- 1 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- rgamma(n, shape = 3, rate = 3 / exp(eta)); d$w15 <- 1.5
r <- suppressWarnings(frm(y | weights(w15) ~ x + (1 | gr(obs, dist = "student")),
         family = Gamma(link = "log"), data = d))
cat("base: logLik", format(as.numeric(logLik(r)), digits = 10), "shape",
    format(exp(r$estimates$betad[1]), digits = 4), "\n")
