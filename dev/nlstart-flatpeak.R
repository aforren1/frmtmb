# The flat-peak fixture of test-diagnostics-ux.R and test-se-flags.R:
# logLik with the sign probe off and on. Seed 5 as in the fixture.
source("dev/nlstart-common.R")
set.seed(5)
w <- seq(1, 45, by = 0.25)
d <- do.call(rbind, lapply(seq_len(10), function(i) {
  mu_i <- stats::rnorm(1, 10, 1.1)
  logf <- stats::rnorm(1, 1.2, 0.3) - 1.5 * log(w) +
    1.1 * exp(-0.5 * ((w - mu_i) / 1.6)^2)
  data.frame(id = i, w = w, logw = log(w),
             I = stats::rexp(length(w), rate = 1 / exp(logf)))
}))
d$id <- factor(d$id)
f <- bf(I ~ apo - chi * logw + exp(lamp) * exp(-0.5 * ((w - pk) / exp(lsig))^2),
        apo ~ 1 + (1 | id), chi ~ 1, pk ~ 1, lamp ~ 1, lsig ~ 1, nl = TRUE)
for (probe in c(FALSE, TRUE)) {
  options(frmtmb.nl_sign_probe = probe)
  cat("probe", probe, ":", ll(tt(frm(f, family = exponential(link = "log"),
                                     data = d))), "\n")
}
