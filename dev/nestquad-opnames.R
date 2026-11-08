# Operation names on tiny tapes of a few densities, to see which carry
# atomic operators.
library(RTMB)
y <- c(0.3, 0.5)
fs <- list(
  dt = function(p) sum(dt((y - p[1]) / exp(p[2]), exp(p[2]), log = TRUE)),
  lgamma = function(p) sum(lgamma(exp(p[2]) + y)),
  pbeta = function(p) sum(pbeta(plogis(y + p[1]), exp(p[2]), 2)),
  dnorm = function(p) sum(dnorm(y, p[1], exp(p[2]), log = TRUE)))
for (nm in names(fs)) {
  tp <- MakeTape(fs[[nm]], c(0, 0))
  o <- capture.output(tp$print())
  cat(nm, ":", names(table(sub("^[[:space:]]*([^[:space:]]+).*", "\\1",
                                  o[-1]))), "\n")
}
