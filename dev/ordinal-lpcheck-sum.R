# Summarize dev/ordinal-log-lpcheck.txt into the numbers the findings
# quote: rows run, identities held, the largest relative constant and
# the largest gradient. Output: dev/ordinal-log-lpcheck-sum.txt
x <- readLines("dev/ordinal-log-lpcheck.txt")
rows <- sub("^== (.*) ==$", "\\1", grep("^== .* ==$", x, value = TRUE))
rows <- setdiff(rows, grep("^SUMMARY", rows, value = TRUE))
lp <- grep("^LPCHECK", x, value = TRUE)
num <- function(key) {
  as.numeric(sub(paste0("^.*", key, " ([^ ]+).*$"), "\\1", lp))
}
const <- num("const")
grad <- num("grad")
ours <- num("ours")
err <- grep("^ERROR", x)
cat("rows run:", length(rows), "\n")
cat("identities held (LPCHECK lines):", length(lp), "\n")
cat("rows that stopped with an error:", length(err), "\n")
cat("largest |const| / max(1, |logLik|):",
    format(max(abs(const) / pmax(1, abs(ours))), digits = 3), "\n")
cat("largest brms gradient at frmtmb's optimum:",
    format(max(grad), digits = 3), "\n")
