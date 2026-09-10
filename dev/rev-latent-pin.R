# Reviewer, punch round 2: is the shrink weight pinned by anything?
#
# The weight was this round's blocker and it moved twice, 0.5 to 0.9.
# `grep` finds it in one place, `R/lca.R:304` as
# `p <- 0.1 * p + 0.9 * pool`, and nothing in `test-lca.R` mentions the
# shrink at all. So an edit of that constant, in either direction, would
# leave the suite green.
#
# A pin does not need a 200-replicate sweep. Seed 20270476 discriminates
# on its own: this review measured w = 0, 0.25, 0.5 and 0.99 all losing
# it by 242 to 256 log-likelihood units while 0.75, 0.9 and 0.95 reach
# poLCA's optimum. This checks that the SHIPPED default reaches it, and
# what an assertion on that seed would cost the ordinary suite.
#
#   Rscript dev/rev-latent-pin.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env3.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
source(paste0(D, "latent-lca-sim.R"))
suppressMessages(loadNamespace("poLCA"))
K <- 4L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))
SEED <- 20270476L
POLCA <- -11507.6035576          # adjudicated in dev/rev-latent-oos3.R

s <- lca_sim(seed = SEED)
t0 <- Sys.time()
fit <- suppressWarnings(suppressMessages(
  frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
      control = tight)))
sf <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
ll <- as.numeric(logLik(fit))
cat("seed", SEED, "\n")
cat("shipped default logLik:", format(ll, digits = 12),
    sprintf("  (%.2f s)\n", sf))
cat("poLCA(nrep = 10)      :", format(POLCA, digits = 12), "\n")
cat("relative gap          :", format(abs(ll - POLCA) / abs(ll),
                                      digits = 4), "\n")
cat("reaches it            :", abs(ll - POLCA) < 1e-6 * abs(ll), "\n")
cat("\nwhat this seed separates, from dev/rev-latent-oos4.R:\n")
cat("  w = 0, 0.25, 0.5, 0.99 lose it by 242 to 256 units\n")
cat("  w = 0.75, 0.9, 0.95 reach it\n")
cat("so one lca() fit on this seed, at", format(sf, digits = 3),
    "s, is a pin on the constant\n")
cat("in BOTH directions. Nothing in test-lca.R does that today.\n")
cat("\ngrep for the constant in the package:\n")
src <- readLines(paste0(
  "C:/Users/adf44/source/r/frmtmb-wt-latent/extensions/frmtmb.latent",
  "/R/lca.R"))
hits <- grep("0[.]9 [*] pool|0[.]1 [*] p", src)
for (h in hits) cat(sprintf("  R/lca.R:%d  %s\n", h, trimws(src[h])))
tst <- readLines(paste0(
  "C:/Users/adf44/source/r/frmtmb-wt-latent/extensions/frmtmb.latent",
  "/tests/testthat/test-lca.R"))
cat("mentions of 'shrink' or 'pool' in test-lca.R:",
    length(grep("shrink|pool", tst)), "\n")
