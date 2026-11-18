# Per-parameter MLE against the stanli posterior, from dev/samplers-q1.R's
# records; the MLE itself is recovered as median + z * sd.
ids <- commandArgs(trailingOnly = TRUE)
for (id in ids) {
  r <- readRDS(file.path("dev/samplers-q1-out", paste0(id, ".rds")))
  s <- r$sum_stanli
  z <- r$z_mle
  h <- names(z)[!grepl("^(r_|z_|zgp_|s_|lp__|lprior)", names(z))]
  m <- s$median[match(h, s$variable)]
  sd <- s$sd[match(h, s$variable)]
  cat("==", id, r$key, "\n")
  print(data.frame(variable = h, mle = signif(m + z[h] * sd, 4),
                   post_median = signif(m, 4), post_sd = signif(sd, 3),
                   z = round(z[h], 2)), row.names = FALSE)
  if (length(r$mle_missing)) cat("MLE names not in the posterior:",
                                 head(r$mle_missing, 8), "\n")
}
