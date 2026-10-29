# Lane meme: z-scores of frm_sample() against brms 2.23.0 from the
# summaries dev/meme-agree.R writes. Output: dev/meme-out/agree-cmp.txt
out <- "C:/Users/adf44/source/r/frmtmb-wt-meme/dev/meme-out"
b <- readRDS(file.path(out, "agree-brms.rds"))
for (arm in c("lane", "base")) {
  f <- file.path(out, paste0("agree-", arm, ".rds"))
  if (!file.exists(f)) next
  s <- readRDS(f)
  v <- intersect(rownames(b), rownames(s))
  z <- cbind(
    mean = (s[v, "mean"] - b[v, "mean"]) /
      sqrt(s[v, "mcse"]^2 + b[v, "mcse"]^2),
    sd = (s[v, "sd"] - b[v, "sd"]) /
      sqrt(s[v, "mcse_sd"]^2 + b[v, "mcse_sd"]^2),
    q05 = (s[v, "q05"] - b[v, "q05"]) /
      sqrt(s[v, "mcse_q05"]^2 + b[v, "mcse_q05"]^2),
    q95 = (s[v, "q95"] - b[v, "q95"]) /
      sqrt(s[v, "mcse_q95"]^2 + b[v, "mcse_q95"]^2))
  cat("==== arm", arm, ": z = (frmtmb - brms) / sqrt(mcse^2 + mcse^2)\n")
  print(round(z, 2))
  cat("missing in arm:", setdiff(rownames(b), rownames(s)), "\n")
  cat(sprintf("max |z| %.2f over %d variables x 4 statistics\n",
              max(abs(z)), length(v)))
}
