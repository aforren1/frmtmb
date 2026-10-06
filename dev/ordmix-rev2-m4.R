# Reviewer of lane ordmix, re-check of m4: the lane found the taped
# gradient of ord_log_interior() with and without the gap floor differ
# on 8675 of 200000 entries, up to 4e-3 relative (dev/ordmix-p1-m4.R).
# Which form is right? Both against the analytic derivative of
# f(a, b) = log(F(a) - F(b)), F the logistic, written stably:
#   df/da =  F'(a) / (F(a) - F(b)),  df/db = -F'(b) / (F(a) - F(b)),
#   log F'(x) = log F(x) + log F(-x),
#   log(F(a) - F(b)) = -b + log(-expm1(b - a)) + log F(a) + log F(b).
# Same pairs as the lane's script (seed 20261005).
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
oli <- frmtmb:::ord_log_interior
set.seed(20261005)
n <- 100000
b <- runif(n, -40, 40)
a <- b + exp(runif(n, log(1e-12), log(50)))
x <- c(a, b)
gsum <- function(on) {
  RTMB::MakeTape(function(p) sum(oli(p[seq_len(n)], p[n + seq_len(n)], on)),
                 x)
}
g_on <- as.numeric(gsum(TRUE)$jacobian(x))
g_off <- as.numeric(gsum(FALSE)$jacobian(x))
lF <- function(z) plogis(z, log.p = TRUE)
ld <- -b + log(-expm1(b - a)) + lF(a) + lF(b)
ga <- exp(lF(a) + lF(-a) - ld)
gb <- -exp(lF(b) + lF(-b) - ld)
g_ref <- c(ga, gb)
rel <- function(g) abs(g - g_ref) / abs(g_ref)
d <- which(g_on != g_off)
cat(sprintf("entries where the two tapes differ: %d of %d\n", length(d),
            2 * n))
for (nm in c("on", "off")) {
  g <- if (nm == "on") g_on else g_off
  r <- rel(g)
  cat(sprintf("floor %-3s vs analytic, all entries: median rel %.3g, 99%% %.3g, max %.3g\n",
              nm, median(r), quantile(r, 0.99), max(r)))
  cat(sprintf("floor %-3s vs analytic, where the tapes differ: median %.3g, 99%% %.3g, max %.3g\n",
              nm, median(r[d]), quantile(r[d], 0.99), max(r[d])))
}
better <- sum(rel(g_on)[d] < rel(g_off)[d])
worse <- sum(rel(g_on)[d] > rel(g_off)[d])
cat(sprintf("where they differ, the floored tape is closer in %d, farther in %d, tied %d\n",
            better, worse, length(d) - better - worse))
big <- d[pmax(rel(g_on)[d], rel(g_off)[d]) > 1e-6]
cat(sprintf("entries with either error > 1e-6: %d; which side: b in %d; gap range %.3g..%.3g; b range %.3g..%.3g\n",
            length(big), sum(big > n),
            min((a - b)[(big - 1) %% n + 1]), max((a - b)[(big - 1) %% n + 1]),
            min(b[(big - 1) %% n + 1]), max(b[(big - 1) %% n + 1])))
k <- big[order(-pmax(rel(g_on)[big], rel(g_off)[big]))][1:6]
i <- (k - 1) %% n + 1
print(data.frame(a = a[i], b = b[i], wrt = ifelse(k > n, "b", "a"),
                 analytic = g_ref[k], floor_on = g_on[k], floor_off = g_off[k]))
# the size of these derivatives against the other partial of the pair:
# an error that matters for an optimizer is relative to the gradient
# vector, not to a 1e-17 entry
j <- (big - 1) %% n + 1
scale <- pmax(abs(ga[j]), abs(gb[j]))
cat(sprintf("largest error relative to the pair's larger partial: on %.3g, off %.3g\n",
            max(abs(g_on[big] - g_ref[big]) / scale),
            max(abs(g_off[big] - g_ref[big]) / scale)))
