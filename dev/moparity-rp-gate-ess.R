# Which statistics have the smallest ESS in dev/moparity-rp-gate.R output.
#   Rscript dev/moparity-rp-gate-ess.R <glob>
fs <- Sys.glob(commandArgs(TRUE)[1])
z <- do.call(rbind, lapply(fs, function(f) cbind(readRDS(f), file = basename(f))))
o <- z[order(pmin(z$ess_c, z$ess_n)), ]
print(head(o[, c("file", "seed", "var", "z_mean", "z_sd", "ess_c", "ess_n")], 12),
      digits = 3)
