# Summary of dev/famuni-r2-xiedge.R's two runs: every fit inside the
# 1e-3 band, with its slope and what gev_fit_check() said, and whether
# any fit outside the band warned or messaged (none should).
# Output: dev/famuni-out/r2-xiedge-sum.txt
for (l in c("log1p", "identity")) {
  r <- readRDS(file.path("C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-out",
                         paste0("r2-xiedge-", l, ".rds")))
  band <- r$xi < -1 + 1e-3
  steep <- band & abs(r$dll) > 0.1
  flat <- band & abs(r$dll) <= 0.1
  cat(sprintf(paste0("\n== %s: %d fits; in the band %d; warned %d (all in ",
                     "band: %s); messaged %d (all in band: %s)\n"),
              l, nrow(r), sum(band), sum(r$warn), all(band[r$warn]),
              sum(r$msg), all(band[r$msg])))
  cat(sprintf(paste0("   in band, |d logLik / d xi| > 0.1: %d, all warned: ",
                     "%s; <= 0.1: %d, all messaged: %s\n"),
              sum(steep), all(r$warn[steep]), sum(flat), all(r$msg[flat])))
  cat("   warned fits, d logLik / d xi range:",
      format(range(r$dll[r$warn]), digits = 3), "\n")
  print(r[band & (flat | r$dll > 0),
          c("xi0", "n", "seed", "xi", "dll", "conv", "warn", "msg")],
        row.names = FALSE)
}
