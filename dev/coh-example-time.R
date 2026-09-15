## How long the `?cross_wishart` example takes, in its own process.
##
## R CMD check reported a NOTE for it at 4.69 + 0.59 seconds of CPU
## against a 5 second threshold. This lane did not touch the example
## (the Rd diff is insertions only, none inside \examples), so the
## number is worth having beside the claim rather than the claim alone.
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
t <- system.time({
  set.seed(4)
  src <- rnorm(4096)
  xs <- frm_cross_spectrum(src + rnorm(4096), 0.9 * src + rnorm(4096),
                           segments = 16)
  xs <- xs[xs$freq < 0.1, ]
  fit <- frmtmb::frm(
    frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
               pow2 ~ 1, coh ~ 1, phase ~ 1),
    family = cross_wishart(), data = xs)
  print(frm_coherence(fit)[1, ])
})
print(t)
cat("user + system =", t[["user.self"]] + t[["sys.self"]], "\n")
