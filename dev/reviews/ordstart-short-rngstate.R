# Short check of lane ordstart, finding 1: mixture_starts_draws() puts
# back the caller's RNG kinds and .Random.seed, and leaves a missing seed
# missing.
#   Rscript ordstart-short-rngstate.R
.libPaths(c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
dr <- frmtmb:::mixture_starts_draws
ge <- globalenv()
if (exists(".Random.seed", envir = ge)) rm(".Random.seed", envir = ge)
z0 <- dr(3, 1)
cat("missing seed stays missing:", !exists(".Random.seed", envir = ge), "\n")
for (k in list(c("Mersenne-Twister", "Inversion", "Rejection"),
               c("L'Ecuyer-CMRG", "Box-Muller", "Rejection"),
               c("Knuth-TAOCP-2002", "Ahrens-Dieter", "Rounding"))) {
  suppressWarnings(RNGkind(k[1], k[2], k[3]))
  set.seed(5)
  before <- get(".Random.seed", envir = ge)
  kb <- RNGkind()
  z <- dr(3, 1)
  cat(sprintf("kind %-16s draws same %s | kinds back %s | seed back %s\n",
              k[1], identical(z, z0), identical(RNGkind(), kb),
              identical(get(".Random.seed", envir = ge), before)))
}
