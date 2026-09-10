# Lane `latent`, punch round 2: the IN-SAMPLE grid extended to the
# right-hand end.
#
# `dev/latent-2p4-shrink.tsv` holds w = 0, 0.1, 0.25, 0.5, 0.75, 0.9 on
# the lane's own 200 seeds. Locating the plateau's upper edge needs
# 0.95, 0.99 and 1 as well, on those same seeds, so that the in-sample
# and out-of-sample grids can be read side by side without pooling.
#
# The poLCA reference is READ from `dev/latent-2p4-lca.tsv`. No EM runs
# here.
#
#   Rscript dev/latent-2p4-shrink2.R [nseed]
#
# Seeds: the lane's own, 20260910 + 0..199.

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("dev/latent-lca-sim.R")
source("dev/latent-2p4-startfix-fns.R")

args <- commandArgs(trailingOnly = TRUE)
NSEED <- if (length(args)) as.integer(args[[1L]]) else 200L
WS <- c(0.95, 0.99, 1)
K <- 4L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))
ref <- utils::read.delim("dev/latent-2p4-lca.tsv", check.names = FALSE)
OUT <- "dev/latent-2p4-shrink2.tsv"
wnm <- paste0("ll_w", sub("[.]", "", format(WS)))
cat(paste(c("seed", "ll_polca", wnm), collapse = "\t"), "\n",
    sep = "", file = OUT)

for (i in seq_len(min(NSEED, nrow(ref)))) {
  seed <- ref$seed[i]
  s <- lca_sim(seed = seed)
  nc <- rep(2L, ncol(s$Y))
  sl <- slice_kmeans(s$Y, nc, K)
  lls <- vapply(WS, function(w) {
    f <- try(suppressWarnings(suppressMessages(
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
          control = tight,
          start = profiles_from(s$Y, nc, K, sl, w = w)))), silent = TRUE)
    if (inherits(f, "try-error")) NA_real_ else as.numeric(logLik(f))
  }, numeric(1))
  cat(paste(c(seed, formatC(c(ref$ll_polca[i], lls), digits = 6,
                            format = "f")),
            collapse = "\t"), "\n", sep = "", file = OUT, append = TRUE)
  if (i %% 25L == 0L) { cat("  ", i, " seeds\n", sep = ""); flush(stdout()) }
}
cat("wrote", OUT, "\n")
