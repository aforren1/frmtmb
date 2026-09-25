# Merge the six chunks of dev/phase3b-cdf-reference-v3.R and report how
# far the two independent 700-bit routes agree.
# Output: dev/phase3b-log/cdf-reference-v3.csv and cdf-reference-v3.txt,
# and the test fixture extensions/frmtmb.eam/tests/testthat/fixtures/
# wiener-rtcdf-ref.csv.
x <- do.call(rbind, lapply(1:6, function(k) {
  utils::read.csv(sprintf("dev/phase3b-log/cdf-ref3-chunk%d.csv", k),
                  colClasses = c(lS = "character", lF = "character",
                                 lFl = "character", lFu = "character"))
}))
x <- x[order(x$v, x$a, x$w, x$u), ]
stopifnot(nrow(x) == 3150L, !anyDuplicated(x[, c("v", "a", "w", "u")]))
for (nm in c("lS", "lF", "lFl", "lFu")) {
  stopifnot(!anyNA(as.numeric(x[[nm]])))
}
utils::write.csv(x, "dev/phase3b-log/cdf-reference-v3.csv", row.names = FALSE)
utils::write.csv(x[, c("v", "a", "w", "u", "t", "lS", "lF", "lFl", "lFu")],
                 "extensions/frmtmb.eam/tests/testthat/fixtures/wiener-rtcdf-ref.csv",
                 row.names = FALSE)
q <- function(z) format(max(z, na.rm = TRUE), digits = 3)
out <- c(sprintf("700-bit reference, %d rows; |v| a up to %g", nrow(x),
                 max(abs(x$v) * x$a)),
         sprintf("identity |Fl + Fu + S - 1| on the images, max: %s", q(x$ident)),
         sprintf("images against eigenfunctions, u >= 0.05 (%d rows), max relative disagreement: S %s, Fl %s, Fu %s",
                 sum(!is.na(x$agreeS)), q(x$agreeS), q(x$agreeFl), q(x$agreeFu)),
         sprintf("log S from %.1f to %.3g", min(as.numeric(x$lS)),
                 max(as.numeric(x$lS))),
         sprintf("log F_lower from %.1f; log F_upper from %.1f",
                 min(as.numeric(x$lFl)), min(as.numeric(x$lFu))))
writeLines(out, "dev/phase3b-log/cdf-reference-v3.txt")
cat(out, sep = "\n")
