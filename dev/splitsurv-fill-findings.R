# Paste generated blocks into dev/splitsurv-findings.md at their markers
# (<!-- splitsurv:<name> -->), so that no number in it is typed by hand.
#   Rscript dev/splitsurv-fill-findings.R
f <- "dev/splitsurv-findings.md"
x <- readLines(f)
blocks <- c(compat = "dev/splitsurv-out/compare.txt",
            curve = "dev/splitsurv-out/compare.txt",
            fits = "dev/splitsurv-out/fit-compare.txt",
            grep = "dev/splitsurv-out/grep.txt",
            tests = "dev/splitsurv-out/tests.md",
            checks = "dev/splitsurv-out/checks.txt")
pick <- function(nm, lines) {
  if (nm == "curve") {
    lines[seq(grep("^== curve outputs", lines), grep("^== frm_compat", lines) - 1L)]
  } else if (nm == "compat") {
    lines[grep("^== frm_compat", lines):length(lines)]
  } else lines
}
for (nm in names(blocks)) {
  mk <- paste0("<!-- splitsurv:", nm, " -->")
  i <- which(x == mk)
  if (length(i) != 1L || !file.exists(blocks[[nm]])) next
  body <- pick(nm, readLines(blocks[[nm]]))
  # a block already pasted runs to the closing marker
  end <- which(x == paste0("<!-- /splitsurv:", nm, " -->"))
  rest <- if (length(end)) x[(end + 1L):length(x)] else x[(i + 1L):length(x)]
  fence <- if (nm == "tests") character(0) else "```"
  x <- c(x[1:i], fence, body, fence, paste0("<!-- /splitsurv:", nm, " -->"),
         rest)
}
writeLines(x, f)
