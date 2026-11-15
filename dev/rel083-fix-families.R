# Merge step 1 (famend x silent2), R/families.R
# mixture_ord_degenerate_check(): the three conflicts resolved as silent2's
# findings ("Merge note for consolidation (famend, R/families.R)") say:
# keep both sides, "flagged in five cases", silent2's flags and famend's
# `shut`, `bad <- step | runoff | unplaced | dg$collapsed %in% TRUE | shut`,
# silent2's collapsed text, then famend's shut branch.
f <- "R/families.R"
x <- readLines(f)
hunks <- function(x) {
  a <- grep("^<<<<<<< rel$", x); m <- grep("^=======$", x)
  b <- grep("^>>>>>>> silent2$", x)
  stopifnot(length(a) == 3, length(m) == 3, length(b) == 3)
  list(a = a, m = m, b = b)
}
h <- hunks(x)
side <- function(k, s) {
  if (s == "rel") x[(h$a[k] + 1):(h$m[k] - 1)] else x[(h$m[k] + 1):(h$b[k] - 1)]
}
# 1: silent2's two bullets, then famend's cs() bullet (the rel side less
# its first, collapsed, bullet)
r1 <- side(1, "rel"); s1 <- side(1, "s2")
cs_at <- grep("^#' - A `cs[(][)]` effect has closed", r1)
stopifnot(length(cs_at) == 1)
res1 <- c(s1, r1[cs_at:length(r1)])
# 2: silent2's flags, famend's shut, the one bad
r2 <- side(2, "rel"); s2 <- side(2, "s2")
stopifnot(grepl("^  shut <- ", r2[1]), grepl("^  bad <- ", s2[length(s2)]))
res2 <- c(s2[-length(s2)], r2[1],
          "  bad <- step | runoff | unplaced | dg$collapsed %in% TRUE | shut")
# 3: silent2's collapsed text, then famend's shut branch
r3 <- side(3, "rel"); s3 <- side(3, "s2")
sh <- grep("^      if [(]shut[[]i[]][)] [{]$", r3)
stopifnot(length(sh) == 1)
res3 <- c(s3, r3[(sh - 1):length(r3)])
out <- c(x[1:(h$a[1] - 1)], res1, x[(h$b[1] + 1):(h$a[2] - 1)], res2,
         x[(h$b[2] + 1):(h$a[3] - 1)], res3, x[(h$b[3] + 1):length(x)])
fc <- grep("flagged in four cases", out)
cat("four-cases lines:", fc, "\n")
stopifnot(length(fc) == 1)
out[fc] <- sub("flagged in four cases", "flagged in five cases", out[fc])
con <- file(f, "wb"); writeLines(out, con, sep = "\r\n"); close(con)
cat("resolved\n")
