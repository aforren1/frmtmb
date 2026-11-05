# Fill the short check of fix round 2 from the logs and append it to the review once.
S <- "C:/Users/adf44/AppData/Local/Temp/2/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/r2-section.md"
s <- paste(readLines(S, warn = FALSE), collapse = "\n")
O <- "dev/review-predread-out/r2/"
st <- vapply(c("frmtmb", "frmtmb.sample", "frmtmb.latent"), function(p) {
  x <- readLines(paste0(O, "rcheck/", p, "-check.txt"), warn = FALSE)
  sub("^Status: ", "", tail(grep("^Status:", x, value = TRUE), 1L))
}, "")
stopifnot(all(st == "1 NOTE"))
sub1 <- function(s, old, new) {
  parts <- strsplit(s, old, fixed = TRUE)[[1L]]
  if (length(parts) != 2L) stop("matched ", length(parts) - 1L, ": ", old)
  paste0(parts[1L], new, parts[2L])
}
s <- sub1(s, "CHECKSTATUS", "core, frmtmb.sample and frmtmb.latent each `Status: 1 NOTE`, the V8 math-rendering one.")
s <- sub1(s, "No test and no audit cell holds more than about 400 new rows, so 474 green files say nothing about it. The merged tree carries it too.",
          "The suite is green because no test reads these on that many new rows, so 474 files say nothing about it. The merged tree carries it too (`rtests-merge.log`).")
s <- sub1(s, "takes 0.9 s on poisson, 10 s on categorical and 41 s on cumulative", "takes 1.1 s on poisson, 10 s on categorical and 41 s on cumulative")
s <- sub1(s, "168 fits, on bernoulli, gaussian, lognormal and poisson, with 8 or 30 groups, 10 or 40 rows per group, group SD 0.2 or 0.8, `(1 | g)` or `(1 + x | g)`.",
          "168 fits, on bernoulli, gaussian, lognormal and poisson, by 14 designs and 3 seeds: 8 or 30 groups, 10 or 40 rows per group, group SD 0.2 or 0.8, `(1 | g)` or `(1 + x | g)`, without 8 groups of 40 at SD 0.2.")
s <- sub1(s, "`avg_comparisons(factor)`", "`avg_comparisons()`, a factor")
s <- sub1(s, "`avg_comparisons(z)`", "`avg_comparisons()`, a second covariate")
s <- sub1(s, "I rule it a release fix, because this round made them worse than the plain Wald it replaced, and the fix is small.",
          "I rule it a release fix, because this round moved them away from brms and the fix is small: a plain Wald at these moments reads 0.068 to 0.275 and 0.143 to 0.293, where the lane reads 0.096 to 0.312 and 0.171 to 0.333.")
s <- sub1(s, "Categorical against frmtmb's draws, category 2: lane 0.230 to 0.397, F1 0.181 to 0.323, draws 0.221 to 0.308.",
          "Categorical against frmtmb's draws, category 2: lane 0.230 to 0.397, F1 0.181 to 0.323, draws 0.221 to 0.308. Its lower end is still wide, because the hump's coefficient share is large on the logit, and it errs on the wide side.")
s <- sub1(s, "a multinomial of 20 trials reads inside the fixed-coefficient draws: 3.67 to 9.58 against 3.60 to 9.68.",
          "a multinomial of 20 trials reads 3.67 to 9.58 against 3.60 to 9.68 for draws at fixed coefficients.")
s <- sub1(s, "base 20 of 20, 14 of 14 caught, with the same differ counts as round 1, cell for cell.",
          "base 20 of 20, 14 of 14 caught, each with the same differ count as round 1.")
s <- sub1(s, "one measured verdict change, p013 `cs_pred()` x `frm_curve` from defect to equivalent, no change outside the lane's set, headline 83.6.",
          "one measured verdict change, p013 `cs_pred()` x `frm_curve` from defect to equivalent, no level changes outside the lane's set, only `seconds` moved among the numeric cells, headline 83.6.")
s <- sub1(s, "All merged readers stop in the same way the lane's do on F2 and F4 until those land.",
          "The merge carries F2's and F4's defects as the lane does, and the four tests fail on it too (`na-merge.log`, `rtests-merge.log`) and pass with the fixes (`rtests-merge-final.log`).")
s <- sub1(s, "A list for the memo. A fit reads about 130 nodes, and `identical()` leaves at the first value that differs.",
          "A list for the memo. A fit reads a few hundred nodes at most, and `identical()` leaves at the first value that differs.")
s <- sub1(s, "bytes, which about 500 rows pass. A fit reads about 130 nodes", "bytes, which about 500 rows pass. A fit reads a few hundred nodes at most")
stopifnot(!grepl("\u2014", s), !grepl("CHECKSTATUS", s))
f <- "dev/reviews/2026-11-02-predread.md"
x <- readLines(f, warn = FALSE)
if (any(x == "## Short check of fix round 2")) stop("already appended")
writeLines(c(x, strsplit(s, "\n", fixed = TRUE)[[1L]]), f)
cat("appended\n")
