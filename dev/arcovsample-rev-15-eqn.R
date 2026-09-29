# REVIEW re-check, detail: where the braces and backslashes in the
# rendered text actually are (so the count in script 14 is not read as
# a broken \eqn), whether \eqn survives Rd2HTML and Rd2latex too, and
# the exact nu at which the OLD Student-t head stops being a number.
#
#   Rscript dev/arcovsample-rev-15-eqn.R

WT <- "C:/Users/adf44/source/r/frmtmb-wt-arcovsample"
files <- c("extensions/frmtmb.sample/man/sample-log_lik.Rd",
           "man/frmtmb-autocor.Rd")

cat("=== source \\eqn uses in the changed Rd\n")
for (f in files) {
  src <- readLines(file.path(WT, f), warn = FALSE)
  k <- grep("eqn", src, fixed = TRUE)
  cat("---- ", basename(f), ": ", length(k), " line(s) with \\eqn\n",
      sep = "")
  for (i in k) cat("  ", src[i], "\n", sep = "")
}

cat("\n=== every rendered line holding a brace or a backslash\n")
for (f in files) {
  txt <- capture.output(tools::Rd2txt(file.path(WT, f), out = stdout(),
                                      package = "x"))
  k <- grep("[\\\\{}]", txt)
  cat("---- ", basename(f), ": ", length(k), "\n", sep = "")
  for (i in k) cat("  [", i, "] ", txt[i], "\n", sep = "")
}

cat("\n=== the same topics through Rd2HTML and Rd2latex\n")
for (f in files) {
  for (fn in c("Rd2HTML", "Rd2latex")) {
    o <- capture.output(do.call(getExportedValue("tools", fn),
                                list(file.path(WT, f), out = stdout(),
                                     package = "x")))
    one <- paste(o, collapse = "\n")
    hit <- grep("max(p, q) + 1", o, fixed = TRUE)
    hit2 <- grep("\\max", o, fixed = TRUE)
    cat("  ", basename(f), " ", fn, ": lines ", length(o),
        "; text form 'max(p, q) + 1' on ", length(hit),
        " line(s); latex form '\\max' on ", length(hit2), " line(s)\n",
        sep = "")
    if (fn == "Rd2latex" && length(hit2)) {
      cat("    ", trimws(o[hit2[1L]]), "\n", sep = "")
    }
    if (fn == "Rd2HTML" && length(hit)) {
      cat("    ", trimws(o[hit[1L]]), "\n", sep = "")
    }
  }
}

cat("\n=== exact nu where the as-written Student-t head stops being a number\n")
old_head <- function(nu, K) lgamma((nu + K) / 2) - lgamma(nu / 2)
for (K in 2:3) {
  lo <- 1e305; hi <- 1e306
  for (it in 1:200) {
    mid <- sqrt(lo * hi)
    if (is.finite(old_head(mid, K))) lo <- mid else hi <- mid
  }
  cat("  K = ", K, ": finite up to nu = ", format(lo, digits = 6),
      ", NaN from nu = ", format(hi, digits = 6), "\n", sep = "")
}
cat("  lgamma() overflows to Inf above argument ",
    format(local({ lo <- 1e305; hi <- 1e306
      for (it in 1:200) { mid <- sqrt(lo * hi)
        if (is.finite(lgamma(mid))) lo <- mid else hi <- mid }
      lo }), digits = 6), "\n", sep = "")
cat("  the NEWS entry names 3.6e305\n")
cat("  old head at 3.6e305, K = 2: ", format(old_head(3.6e305, 2)),
    " (the truth is about ", format((2 / 2) * log(3.6e305 / 2),
                                    digits = 6), ")\n", sep = "")
cat("\nDONE\n")
