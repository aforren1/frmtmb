# Reviewer of lane sampleprior: a mixture with a constant dpar
# (sigma1 = 1, a mapped betad entry): the weights' dirichlet(1) must
# still land on the weight's coordinate, not the one the map shifted.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
lt <- function(v, df, mu, s) stats::dt((v - mu) / s, df, log = TRUE) - log(s)
mx <- mean(x)
d <- D$mix
sc <- strsplit(as.character(brms::stancode(
  brms::bf(y ~ x, sigma1 = 1), data = d,
  family = brms::mixture(gaussian, gaussian))), "\n")[[1]]
cat(grep("lprior \\+=", sc, value = TRUE), sep = "\n")
ff <- bf(y ~ x, sigma1 = 1)
fam <- frmtmb::mixture(gaussian(), gaussian())
od <- cap_obj(ff, d, fam)
of <- cap_obj(ff, d, fam, prior = "flat")
print(od$env$par)
set.seed(4)
P <- lapply(1:6, function(k) od$env$par + rnorm(length(od$env$par), 0, 0.4))
fl <- vapply(P, function(q) logprior_at(od, of, q), 1)
bl <- vapply(P, function(q) {
  lt(q[1] + q[2] * mx, 3, 5.3, 2.5) + lt(q[3] + q[4] * mx, 3, 5.3, 2.5) +
    lt(exp(q[5]), 3, 0, 2.5) + q[5] +
    log(plogis(q[6]) * (1 - plogis(q[6])))
}, 1)
cat("max |d frm - d brms|", max(abs((fl - fl[1]) - (bl - bl[1]))), "\n")
