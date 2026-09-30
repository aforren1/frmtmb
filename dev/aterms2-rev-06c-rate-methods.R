# Reviewer, claim 4: post-fit paths that read the mean without going
# through fitted(): residuals(type = "response"), response_mean(),
# emmeans, dharma_residuals, a rate model against its offset(log(time))
# twin at the SAME parameter vector (the rate fit's estimates copied
# into the offset fit). Seed 606 data. Log: dev/aterms2-rev-log-06c.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(606)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4))
d$y <- rpois(n, exp(0.2 + 0.5 * d$x) * d$time)
fr <- q(frm(y | rate(time) ~ x, data = d, family = poisson()))
fo <- q(frm(y ~ x + offset(log(time)), data = d, family = poisson()))
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
cat("residuals response  rel:", rel(residuals(fr)[, 1], residuals(fo)[, 1]), "\n")
cat("residuals pearson   rel:", rel(residuals(fr, type = "pearson")[, 1],
                                    residuals(fo, type = "pearson")[, 1]), "\n")
cat("residuals deviance  rel:", tryCatch(rel(residuals(fr, type = "deviance")[, 1],
                                    residuals(fo, type = "deviance")[, 1]),
                                    error = function(e) conditionMessage(e)), "\n")
rm1 <- tryCatch(response_mean(fr), error = function(e) conditionMessage(e))
rm2 <- tryCatch(response_mean(fo), error = function(e) conditionMessage(e))
cat("response_mean:", if (is.numeric(rm1)) rel(as.numeric(rm1), as.numeric(rm2)) else rm1, "\n")
dh1 <- tryCatch(q(dharma_residuals(fr, n_sim = 200, seed = 1)), error = function(e) conditionMessage(e))
dh2 <- tryCatch(q(dharma_residuals(fo, n_sim = 200, seed = 1)), error = function(e) conditionMessage(e))
if (is.list(dh1) && is.list(dh2)) {
  cat("dharma scaled residuals rel:", rel(dh1$scaledResiduals, dh2$scaledResiduals),
      " fitted rel:", rel(dh1$fittedPredictedResponse, dh2$fittedPredictedResponse), "\n")
} else cat("dharma:", substr(paste(dh1), 1, 150), "|", substr(paste(dh2), 1, 150), "\n")
em1 <- tryCatch(summary(q(emmeans::emmeans(fr, ~ 1, type = "response"))),
                error = function(e) conditionMessage(e))
em2 <- tryCatch(summary(q(emmeans::emmeans(fo, ~ 1, type = "response"))),
                error = function(e) conditionMessage(e))
cat("emmeans rate:\n"); print(em1); cat("emmeans offset:\n"); print(em2)
cat("mean(time)", mean(d$time), " mean(log(time))", mean(log(d$time)), "\n")
