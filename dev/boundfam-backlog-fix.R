# Lane boundfam fix round: move the lane's filed section to just above
# "## Filed at the 0.81.0 release", indent its CLOSED lines, and extend
# the joint-runaway item. Run once:
#   Rscript dev/boundfam-backlog-fix.R
f <- "C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/test-backlog.md"
x <- readLines(f, warn = FALSE)
s <- which(x == "## Filed in the round of 2026-11-11")
stopifnot(length(s) == 1L)
e <- which(startsWith(x, "## Triage, 2026-09-07"))
stopifnot(length(e) == 1L, e > s)
x <- x[-(s:(e - 1L))]
# the blank line before the section was the file's own
sec <- c("## Filed in the round of 2026-11-11",
         "",
         "### From lane boundfam",
         "",
         "- **A joint runaway is still called flat** (dev/boundfam-findings.md,",
         "  \"Not done\"): a poisson fit whose one count is at the largest x",
         "  (`dev/boundfam-repro.R`, case pois1max) runs (Intercept) to -436.5",
         "  and x to 208.8 together, with a code-7 warning that says the stop is",
         "  at the maximum and the \"flat\" warning. glm() warns \"fitted rates",
         "  numerically 0\". The limit rule reads one coefficient at a time, and",
         "  each coefficient alone costs likelihood there. The verdict depends",
         "  on the coding (dev/reviews/2026-11-11-boundfam.md, m2): poisson",
         "  `y ~ g` with group a all zero gives (Intercept) -22.3 and gb 23.5",
         "  \"flat\", while `y ~ 0 + g` gives ga \"limit\" and the warning that the",
         "  maximum likelihood estimate does not exist; the zero-inflated and",
         "  negbinomial forms of `y ~ g` stay \"flat\".",
         "- **A negbinomial shape at its limit stays \"nocurv\" behind a",
         "  convergence warning** (dev/boundfam-findings.md, \"Fix round\"):",
         "  seed 17 of poisson data, log shape 23.2, code 8 and a gradient of",
         "  1.35e-3 above `grad_tol`. The convergence verdict makes the stop",
         "  short, so the standard-error analysis is not applied, and the limit",
         "  rule, which would fire there, is not asked.",
         "- **3022 and 3037 are family fit-end checks, not the boundary rule**",
         "  (dev/boundfam-findings.md, \"Not done\"): both stay open as filed.",
         "")
h <- which(startsWith(x, "## Filed at the 0.81.0 release"))
stopifnot(length(h) == 1L)
x <- c(x[seq_len(h - 1L)], sec, x[h:length(x)])
cl <- startsWith(x, "CLOSED in the round of 2026-11-11 (lane boundfam)")
x[cl] <- paste0("  ", x[cl])
# binary mode: the file has LF endings, and text mode on Windows writes CRLF
con <- file(f, "wb")
writeLines(x, con, sep = "\n")
close(con)
cat("moved; CLOSED lines indented:", sum(cl), "\n")
