# Paste the logs into dev/simnewdata-findings.md verbatim, so every count
# in it is generated rather than typed. Each placeholder line is replaced
# by an indented block of the named log lines.
#   Rscript dev/simnewdata-fill.R
f <- "dev/simnewdata-findings.md"
x <- readLines(f)
L <- function(p) readLines(file.path("dev/simnewdata-log", p), warn = FALSE)
blk <- function(lines) c(paste0("    ", lines))
pick <- function(p, pat) {
  y <- L(p)
  y[grepl(pat, y)]
}
fills <- list(
  "INVARIANT-BLOCK" = blk(substr(L("invariant.txt"), 1, 72)),
  "BRMS-ARMA-BLOCK" = c(
    "`dev/simnewdata-brms-arma.R`, `brms-arma.txt`: brms 2.23.0 and frmtmb",
    "on `y ~ 1 + ar(time, gr = g, cov = TRUE)`, 40 groups of 6, AR(1) 0.7,",
    "newdata rows at times 2 and 3 of one new group and time 2 of another:",
    "", blk(tail(L("brms-arma.txt"), 4))),
  "SEEFAIL-BLOCK" = blk(c(
    pick("seefail-base-test-simulate-newdata.R.txt", "^RESULT"),
    pick("seefail-base-test-pp-check-types.R.txt", "^RESULT"),
    pick("seefail-base-test-draws-methods.R.txt", "^RESULT"),
    pick("seefail-base-test-counterfactual.R.txt", "^RESULT"))),
  "CALLERS-BLOCK" = blk(c(pick("callers-base.txt", "^base:"),
                          pick("callers-lane.txt", "^lane:"),
                          L("callers-compare.txt"))),
  "RPLOTS-BLOCK" = blk(c(L("rplots-base.txt"), L("rplots-lane.txt"))),
  "LEDGER-BLOCK" = c(
    "Before (`before-ledger-summary.md`, the committed 0.62.0 summary) and",
    "after (`dev/brmsport-log/ledger-summary.md`, `ledger-final.txt`):", "",
    blk(c(pick("before-ledger-summary.md", "^Bin 1|^\\| \\*\\*total|brmsfit-methods.R|refuses-accepted"),
          "->",
          pick("ledger-final.txt", "^Bin 1|^\\| \\*\\*total|brmsfit-methods.R` \\||refuses-accepted")))),
  "BITWISE-BLOCK" = blk(tail(L("bitwise-compare.txt"), 2))
)
for (k in names(fills)) {
  i <- which(x == k)
  if (length(i) != 1L) next
  x <- c(x[seq_len(i - 1L)], fills[[k]], x[-seq_len(i)])
}
con <- file(f, "wb")
writeLines(x, con, sep = "\n")
close(con)
cat("left:", grep("-BLOCK$|-SENTENCE$", x, value = TRUE), "\n")
