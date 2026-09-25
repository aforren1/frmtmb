# The recovery arms ran on builds earlier than the final one: `cens` on
# the build with one blend center for the distribution function,
# `contfix` on the build before it. Refit a few of their seeds on the
# FINAL build (dev/phase3b-eam-recovery4.R) and compare.
# Output: dev/phase3b-log/buildcheck.txt
out <- character(0)
for (f in list.files("dev/phase3b-log/recov-final", full.names = TRUE)) {
  b <- readRDS(f)
  a <- readRDS(file.path("dev/phase3b-log/recov", basename(f)))
  if (!is.null(a$error)) {
    a <- readRDS(file.path("dev/phase3b-log/recov-rerun", basename(f)))
  }
  out <- c(out, sprintf(
    "%-14s logLik %.10f -> %.10f (diff %.3g); max |fixef diff| %.3g; identical fixef: %s",
    basename(f), a$loglik, b$loglik, b$loglik - a$loglik,
    max(abs(a$fixef - b$fixef)), identical(a$fixef, b$fixef)))
}
writeLines(out, "dev/phase3b-log/buildcheck.txt")
cat(out, sep = "\n")
