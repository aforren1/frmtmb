## Write the development-version section of NEWS.md, keeping the CRLF
## line endings the file already has and writing no BOM. Idempotent: the
## block is cut back to `# frmtmb 0.64.0` first, so rerunning it after an
## edit replaces the section instead of stacking a second copy.
p <- "C:/Users/adf44/source/r/frmtmb-wt-thresrefit/NEWS.md"
## readLines() does NOT strip the CR here, so it is removed by hand and
## the text-mode connection puts it back on write: adding one by hand
## gave CR CR LF
old <- sub("\r$", "", readLines(p, warn = FALSE))
i <- grep("^# frmtmb 0\\.64\\.0", old)
stopifnot(length(i) >= 1L)
old <- old[i[1L]:length(old)]

new <- c(
"# frmtmb (development version)",
"",
"- A refit inside the package now carries the FITTED model's ordinal",
"  threshold count, and the per-level counts of `thres(gr = )`, instead",
"  of recounting them from the data it is given. `influence()` and",
"  `cooks.distance()` rebuild the design from a subset, so deleting the",
"  last observation in a category used to refit a model with one",
"  threshold fewer: the deleted unit's row of the table held an `NA`, its",
"  Cook's distance was `NA`, and the coefficients above the emptied",
"  category were reported one column early, with nothing but that",
"  trailing `NA` to show it. The deletion now refits the model that was",
"  fitted, and the row agrees with a fit that pins the count by hand.",
"- A refit also carries the fitted model's response CATEGORIES. The model",
"  frame drops a factor level that no row of a subset takes, which",
"  renumbers every category above it, so an ordered-factor response could",
"  give a different influence table from the same data coded as integers.",
"  Emptying an INTERIOR category, and a hand-written `thres(K)`, were the",
"  two cases where it did. Integer, character and ordered-factor codings",
"  now give the same table, bitwise.",
"- `influence()` and `cooks.distance()` now count the deletion refits",
"  that failed. A partly failed table warns with the count and the first",
"  reason; a table whose every refit failed is an error instead of a",
"  silent matrix of `NA`. Two refits cannot represent the fitted model at",
"  all and are refused by name: `data = ` holding a response category, or",
"  a `thres(gr = )` level, that the fit never saw, and `groups = `",
"  deleting a whole `thres(gr = )` level.",
"- `frm_bootstrap()`, `refit()`, `frm_allfit()`, `anova(refit = TRUE)`,",
"  `confint(method = \"profile\")` and the autoscale pre-fit reuse the",
"  assembled design and already kept the count; that is now documented",
"  in `?frm_bootstrap` and pinned by a test. `simulate()` followed by",
"  `frm()`, written out by hand, is not a refit and still counts from the",
"  data it is given, which `?frm_bootstrap` now says.",
"- `frm_simulate(prior = )` used to fail inside `vapply()` with \"values",
"  must be length 1\" whenever one prior covered several parameters at",
"  once. It now says what it cannot do and names the alternative. This is",
"  about the shape of a prior entry, not one family: it is reached by",
"  class `\"Intercept\"` on an ordinal model with more than one threshold,",
"  by class `\"cor\"` on a block with more than one correlation, and by",
"  class `\"ar\"`, `\"ma\"` or `\"cortime\"` above order 1. The underlying",
"  defect, one draw written into every parameter the entry covers, is",
"  refused rather than performed; an entry covering exactly one parameter",
"  still draws.",
"")
writeLines(c(new, old), p)
cat("NEWS.md now starts:\n")
cat(head(readLines(p, warn = FALSE), 3), sep = "\n")
