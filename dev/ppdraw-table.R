# The sweep table of dev/ppdraw-findings.md, generated from the two TSVs
# dev/ppdraw-sweep.R writes. Usage: Rscript dev/ppdraw-table.R
b <- utils::read.delim("dev/ppdraw-sweep-base.tsv", stringsAsFactors = FALSE)
l <- utils::read.delim("dev/ppdraw-sweep-lane.tsv", stringsAsFactors = FALSE)
stopifnot(identical(b$family, l$family), identical(b$case, l$case))
f2 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 2))
cat("| family | fit | z base | z lane | z, no propagation | NA base | NA lane |",
    "held (lane) |\n")
cat("|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(b))) {
  z_b <- if (abs(b$z[i]) > 1e4) formatC(b$z[i], format = "e", digits = 2) else
    f2(b$z[i])
  cat("| ", l$family[i], " | ", l$case[i], " | ", z_b, " | ", f2(l$z[i]),
      " | ", f2(l$z_pe0[i]), " | ", f2(b$na_share[i]), " | ",
      f2(l$na_share[i]), " | ", gsub("`", "", l$held[i]), " |\n", sep = "")
}
big <- function(t) sum(abs(t$z) > 4.5, na.rm = TRUE)
cat("\nfits with |z| > 4.5: base ", big(b), ", lane ", big(l), " of ",
    nrow(l), "\n", sep = "")
cat("fits with a parameter held (lane): ", sum(nzchar(l$held)), "\n", sep = "")
