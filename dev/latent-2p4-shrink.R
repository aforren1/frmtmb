# Lane `latent`, punch round 1: is the shrink weight on a knife edge?
#
# `dev/latent-2p4-startfix.R` measured two weights: a hard k-means start
# (w = 0) fixes the 8 of 200 the shipped start loses and breaks 6 others,
# and shrinking the class profiles halfway to the POOLED profile
# (w = 0.5) reaches poLCA's optimum on 200 of 200. A constant chosen
# from two points is a constant chosen from two points, so this measures
# the whole range and reports where the plateau is.
#
# The reference is READ from dev/latent-2p4-lca.tsv, as before: no new
# poLCA run.
#
#   Rscript dev/latent-2p4-shrink.R [nseed]
#
# Seeds: the sweep's own, 20260910 + 0..199.

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("dev/latent-lca-sim.R")
source("dev/latent-2p4-startfix-fns.R")

args <- commandArgs(trailingOnly = TRUE)
NSEED <- if (length(args)) as.integer(args[[1L]]) else 200L
K <- 4L
WS <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9)
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))
ref <- utils::read.delim("dev/latent-2p4-lca.tsv", check.names = FALSE)
out <- "dev/latent-2p4-shrink.tsv"
cat(paste(c("seed", "ll_polca", paste0("gap_w", sub("[.]", "", WS))),
          collapse = "\t"), "\n", sep = "", file = out)

for (i in seq_len(min(NSEED, nrow(ref)))) {
  seed <- ref$seed[i]
  sm <- lca_sim(seed = seed)
  nc <- rep(2L, ncol(sm$Y))
  sl <- slice_kmeans(sm$Y, nc, K)
  gaps <- vapply(WS, function(w) {
    f <- try(suppressWarnings(suppressMessages(
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = sm$dd,
          control = tight,
          start = profiles_from(sm$Y, nc, K, sl, w = w)))),
      silent = TRUE)
    if (inherits(f, "try-error")) NA_real_ else
      as.numeric(logLik(f)) - ref$ll_polca[i]
  }, numeric(1))
  cat(paste(c(seed, formatC(ref$ll_polca[i], digits = 6, format = "f"),
              formatC(gaps, digits = 6, format = "g")),
            collapse = "\t"), "\n", sep = "", file = out, append = TRUE)
  if (i %% 25L == 0L) { cat("  ", i, " seeds\n", sep = ""); flush(stdout()) }
}

d <- utils::read.delim(out)
tol <- 1e-6 * abs(d$ll_polca)
cat("\nseeds:", nrow(d), "\n")
cat(sprintf("  %-6s %10s %10s\n", "w", "lost", "errored"))
for (k in seq_along(WS)) {
  g <- d[[paste0("gap_w", sub("[.]", "", WS[k]))]]
  cat(sprintf("  %-6s %10d %10d\n", format(WS[k]),
              sum(g < -tol, na.rm = TRUE), sum(is.na(g))))
}
cat("\n(the shipped start loses 8 of", nrow(d), ")\n")
