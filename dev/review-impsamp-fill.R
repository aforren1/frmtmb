# Fill the mutant table into the review and wrap its prose at 80 columns.
p <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/reviews/2026-11-02-impsamp.md"
s <- readLines(p)
tab <- c(
  "`dev/review-impsamp-mutants.R`, `review-impsamp-mutrun.sh`, each mutant a source copy",
  "in `rev-impsamp-src/` and its own library, with the guard that the mutant's library",
  "loaded. An earlier run of mine lacked the guard, installs failed, and every \"mutant\"",
  "passed 277 of 277 on the lane's own build; those logs are not used. `test-importance.R`,",
  "`mut2/test-*.log`, baseline 277 pass:",
  "",
  "| mutant | edit | result |",
  "|---|---|---|",
  "| `partial_rate` | the partial derivative of the curvature, not the total | **277 pass, 0 fail**; one escaped warning |",
  "| `no_cut` | `abs()` in place of the cut of a negative curvature, q = 1 | **277 pass, 0 fail** |",
  "| `half_jacobian` | half the log-determinant's move off the anchor (6 named tests, `mut2/jac-half.txt`) | textbook, curvature, numDeriv, anchor tests all pass; only the gaussian test errors |",
  "| `no_jacobian` | a frozen log-determinant with moving draws (`mut2/jac-no.txt`) | gaussian errors, textbook fails 6 |",
  "| `one_round`, `grad_always` | stop after round 1 whatever the gradient | 276 pass, 1 error: `moves[[2L]]` out of bounds, nothing about accuracy |",
  "| `no_rate` | curvature frozen | gaussian test only, fail 2, error 1 |",
  "| `no_cross` | mode ignores the likelihood's slope | 3 tests fail |",
  "| `no_prior_mode` | the prior's move left out of the mode | 3 tests fail |",
  "| `frozen` | the old estimator | 7 tests fail |",
  "| `map_recycle` | the old `imp_par_list()` | the constant dpar test fails |",
  "| `no_refusal` | the refusal removed | its 2 tests fail |",
  "",
  "`partial_rate` is the lane's third design point, \"the curvature's rate must be the",
  "TOTAL derivative\", with its measurement (4.1 against 8.1, 2 of 648 rising). No test holds",
  "it. `half_jacobian` is a wrong density at every theta but the anchor, and the only test",
  "that sees it is an exact-identity one: nothing compares the objective off the anchor with",
  "an exact value. The whole-file runs of the two Jacobian mutants did not end in 45 minutes,",
  "since every fit climbs the refusal's ladder, and I stopped them; the 6 named tests are what",
  "I have.",
  "",
  "Old tests (`oldtest-importance.txt`, rellib-r16's file on the lane's build): 17 fail, 1",
  "error, in 6 tests, each on the expectation the lane's table names, so the claim that the",
  "fix removed what they constructed holds. See m3 for what the rewrites lost.",
  "",
  "**Fix.** Add two tests. (1) On the poisson slope the lane measured, `obj$he(at)[1, 1]` of",
  "the tape against a central difference of the Laplace gradient, as a ratio, which fails on",
  "`partial_rate`. (2) The objective at theta displaced by 2 SE in each coordinate against",
  "`quadrature = TRUE` or an `integrate()` in the test, within 4 MCSE, which must fail on",
  "`half_jacobian`.")
i <- which(s == "MUTANT_TABLE")
stopifnot(length(i) == 1L)
s <- c(s[seq_len(i - 1L)], tab, s[(i + 1L):length(s)])
j <- grep("| M3 | MUTANTS_M3 | add the tests |", s, fixed = TRUE)
stopifnot(length(j) == 1L)
s[j] <- "| M3 | reverting the total-derivative rate fails no test; half the Jacobian passes all but one | add the tests |"
# wrap prose: a paragraph or a bullet, never a table, heading or indented block
out <- character(0)
k <- 1L
n <- length(s)
while (k <= n) {
  ln <- s[k]
  if (!nzchar(ln) || grepl("^(\\||#|    )", ln)) {
    out <- c(out, ln)
    k <- k + 1L
    next
  }
  e <- k
  while (e < n && nzchar(s[e + 1L]) && !grepl("^(\\||#|- |    )", s[e + 1L])) e <- e + 1L
  txt <- paste(trimws(s[k:e]), collapse = " ")
  bullet <- grepl("^- ", ln)
  w <- strwrap(sub("^- ", "", txt), width = 79, initial = if (bullet) "- " else "",
               prefix = if (bullet) "  " else "")
  out <- c(out, w)
  k <- e + 1L
}
writeLines(out, p)
cat("lines", length(out), " over 80, not tables:",
    sum(nchar(out) > 80 & !grepl("^\\|", out)), "\n")
cat("placeholders left:", sum(grepl("MUTANT|REFUSEC", out)), "\n")
