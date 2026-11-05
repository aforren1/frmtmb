# Review of lane gpgap, item 4: counts and wall clock from the NCP logs.
for (f in Sys.glob("dev/review-gpgap-out/ncp-*.txt")) {
  L <- readLines(f, warn = FALSE)
  L <- L[startsWith(L, "NCP")]
  conv <- as.numeric(sub(".*ncp 1e-12: conv (\\d).*", "\\1", L))
  llc <- as.numeric(sub(".*closed brms ll (\\S+) sdgp.*", "\\1", L))
  lln <- as.numeric(sub(".*ncp 1e-12: conv \\d \\(.*\\) gr \\S+ ll (\\S+) sdgp.*", "\\1", L))
  tc <- as.numeric(sub(".*conv \\d (\\S+)s \\| ncp.*", "\\1", L))
  tn <- as.numeric(sub(".* (\\S+)s \\| chol.*", "\\1", L))
  s <- sub(".* s (\\S+) \\| closed.*", "\\1", L)
  cat(sprintf(paste0("%s: rows %d, relative convergence %d, false convergence %d; ",
                     "at or above the closed form within 0.01: %d; median wall ",
                     "centered %.2f s (with SEs), ncp %.2f s (no SEs); ratio %.1f\n"),
              basename(f), length(L), sum(conv == 0), sum(conv == 1),
              sum(lln >= llc - 0.01), median(tc[-1]), median(tn[-1]),
              median(tn[-1]) / median(tc[-1])))
  tb <- table(s[conv == 1])
  cat("  false convergence by noise:", paste(names(tb), tb, collapse = ", "), "\n")
}
for (f in Sys.glob("dev/review-gpgap-out/ncp2-*.txt")) {
  L <- readLines(f, warn = FALSE)
  for (v in c("rel 1e-8", "rel 1e-10", "rel 1e-12", "abs 1e-12")) {
    x <- grep(paste0("ncp ", v, " "), L, value = TRUE, fixed = TRUE)
    cat(sprintf("%s %-10s conv 0 on %d of %d\n", basename(f), v,
                sum(grepl("conv 0", x, fixed = TRUE)), length(x)))
  }
}
