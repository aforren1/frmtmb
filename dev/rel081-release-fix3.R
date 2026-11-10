# The release review's five minors (dev/reviews/2026-11-10-release.md),
# applied to the release tree. Rerunnable: each edit is applied when its
# old text is there once, skipped when its new text is already there,
# and stops the script otherwise. Line endings are kept as each file has
# them.
#
#   Rscript dev/rel081-release-fix3.R
#
# 1. frmtmb.sample test of the exact reader, through log_lik() with and
#    without the fitted rows as newdata (finding 1).
# 2. loo(), waic() and psis() of draws stop with one classed error on a
#    log_lik() cell that is not finite (finding 2; the coordinator's call:
#    all three stop, as in brms).
# 3. BREAKING marks on code8's bullets 2 to 4 and compat's bullet
#    (finding 3): in dev/rel081-news.R, which this script runs.
# 4. frmtmb.sample NEWS: one floor line, the measured reach of the limit
#    (finding 4): dev/rel081-news.R and dev/test-backlog.md.
# 5. dev/rel081-commit.sh: index sums before the first git call, the
#    staged set checked against take_list, how to resume (finding 5); the
#    brms entries of dev/upstream-bugs.md in numeric order.
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel081"
rd <- function(p) {
  raw <- readChar(p, file.size(p), useBytes = TRUE)
  list(crlf = grepl("\r\n", raw, fixed = TRUE),
       txt = gsub("\r\n", "\n", raw, fixed = TRUE))
}
wrt <- function(f, txt, p) {
  if (f$crlf) txt <- gsub("\n", "\r\n", txt, fixed = TRUE)
  con <- file(p, open = "wb")
  on.exit(close(con))
  writeChar(txt, con, eos = NULL, useBytes = TRUE)
}
n_of <- function(s, pat) lengths(regmatches(s, gregexpr(pat, s, fixed = TRUE)))
edit <- function(path, old, new, what) {
  p <- file.path(REL, path)
  f <- rd(p)
  if (n_of(f$txt, new) >= 1L) {
    cat("already applied:", what, "\n")
    return(invisible(FALSE))
  }
  if (n_of(f$txt, old) != 1L) stop("old text found ", n_of(f$txt, old),
                                   " times: ", what)
  wrt(f, sub(old, new, f$txt, fixed = TRUE), p)
  cat("applied:", what, "\n")
  invisible(TRUE)
}
L <- function(...) paste(c(...), collapse = "\n")

## 2. the classed stop, in frmtmb.sample's loo_matrix() ------------------
loo_r <- "extensions/frmtmb.sample/R/loo.R"
edit(loo_r,
L("#' weights of [psis()], the `loo_*()` functions and `pp_check()`'s",
  "#' `loo_*` types.",
  "#'",
  "#' @param x A `frmtmb_draws` from [frm_sample()], or (for"),
L("#' weights of [psis()], the `loo_*()` functions and `pp_check()`'s",
  "#' `loo_*` types.",
  "#'",
  "#' @section A log-likelihood that is not finite:",
  "#' A cell of the [frmtmb::log_lik()] matrix is `-Inf` when a draw puts",
  "#' an observation outside the support of the response's distribution,",
  "#' for example a `gen_extreme_value()` row below the lower end of a",
  "#' draw's support. A leave-one-out or WAIC estimate does not exist with",
  "#' such a cell. `loo()`, `waic()`, `psis()`, the `loo_*()` functions and",
  "#' `pp_check()`'s `loo_*` types then stop with an error of class",
  "#' `frmtmb_loglik_not_finite` that gives the number of cells, draws and",
  "#' observations. A `NaN` or `+Inf` cell stops them in the same way. brms",
  "#' also stops there: its density is `NaN` outside a support, and the loo",
  "#' package refuses `NaN`. Draws that [frm_sample()] makes put no",
  "#' observation outside its support in practice; a draws object built or",
  "#' edited by hand can.",
  "#'",
  "#' @param x A `frmtmb_draws` from [frm_sample()], or (for"),
"loo.R: help section on a log-likelihood that is not finite")
edit(loo_r,
L("  ll <- draws_as_caller(what, re_formula = FALSE,",
  "                        log_lik(x, ndraws = ndraws, resp = resp, ...))",
  "  # loo::loo.matrix() prints"),
L("  ll <- draws_as_caller(what, re_formula = FALSE,",
  "                        log_lik(x, ndraws = ndraws, resp = resp, ...))",
  "  loo_check_finite(ll, what)",
  "  # loo::loo.matrix() prints"),
"loo.R: loo_matrix() calls loo_check_finite()")
edit(loo_r,
"#' The log-likelihood matrix with the loo package present. `...` reaches",
L("#' Stop on a log-likelihood cell that is not finite, before the loo",
  "#' package sees it. The loo package refuses `+Inf` and `NaN`, and `-Inf`",
  "#' becomes `+Inf` in the importance ratios of `loo()` and `psis()`, but",
  "#' `waic.matrix()` takes `-Inf` and returns an elpd of `NaN` without a",
  "#' word. brms stops in all three, because its density outside a support",
  "#' is `NaN`. One classed error here gives the three functions, and every",
  "#' caller of `loo_matrix()`, the same answer and says which cells.",
  "#'",
  "#' @noRd",
  "loo_check_finite <- function(ll, what) {",
  "  bad <- !is.finite(ll)",
  "  if (!any(bad)) return(invisible(NULL))",
  "  cols <- which(colSums(bad) > 0L)",
  "  rows <- which(rowSums(bad) > 0L)",
  "  unit <- if (is.null(attr(ll, \"unit\"))) \"observations\" else \"units\"",
  "  some <- function(i) {",
  "    paste0(paste(i[seq_len(min(5L, length(i)))], collapse = \", \"),",
  "           if (length(i) > 5L) \", ...\")",
  "  }",
  "  n_minus <- sum(ll[bad] == -Inf, na.rm = TRUE)",
  "  why <- c(",
  "    if (n_minus > 0L) {",
  "      paste0(\"A -Inf cell is an observation that the draw puts outside \",",
  "             \"the support of the response's distribution, where its \",",
  "             \"density is 0.\")",
  "    },",
  "    if (n_minus < sum(bad)) {",
  "      \"A NaN or +Inf cell is a density that is not defined.\"",
  "    })",
  "  frm_stop(what, \": log_lik() is not finite in \", sum(bad), \" cells, at \",",
  "           length(cols), \" of \", ncol(ll), \" \", unit, \" and \", length(rows),",
  "           \" of \", nrow(ll), \" draws (\", unit, \" \", some(cols),",
  "           \"; draws, by row of log_lik(), \", some(rows), \"). \",",
  "           paste(why, collapse = \" \"),",
  "           \" A leave-one-out or WAIC estimate does not exist with such a \",",
  "           \"cell. Leave those draws out with `draw_ids`, or check how the \",",
  "           \"draws were made.\",",
  "           call. = FALSE, class = \"frmtmb_loglik_not_finite\")",
  "}",
  "",
  "#' The log-likelihood matrix with the loo package present. `...` reaches"),
"loo.R: loo_check_finite()")

## 1. the test, appended to test-loo-newdata-fitted.R -------------------
tf <- "extensions/frmtmb.sample/tests/testthat/test-loo-newdata-fitted.R"
f <- rd(file.path(REL, tf))
if (grepl("gev_off_case", f$txt, fixed = TRUE)) {
  cat("already applied: test of the exact reader\n")
} else {
  blk <- rd(file.path(REL, "dev/rel081-fix3-test-block.R"))$txt
  wrt(f, paste0(sub("\n*$", "\n", f$txt), blk), file.path(REL, tf))
  cat("applied: test of the exact reader\n")
}

## 4. the backlog item ---------------------------------------------------
edit("dev/test-backlog.md",
L("- **MINOR. `loo()` and `psis()` of draws stop on a `-Inf` cell.** With",
  "  the exact reader, a draw that puts a `gen_extreme_value()` row outside",
  "  its support gives `-Inf` in `log_lik()`. `loo()` and `psis()` then",
  "  stop with the loo package's \"All input values must be finite or",
  "  -Inf\": the importance log ratio is `-log_lik()`, `+Inf` in that cell,",
  "  and loo 2.10.1's `validate_ll()` refuses `+Inf` (it accepts `-Inf`).",
  "  `waic()` does not stop: its `elpd_waic` is `NaN`. brms 2.23.0's",
  "  `dgen_extreme_value()` gives `NaN` (\"NaNs produced\") in the same",
  "  cell, and `loo()` and `psis()` of that matrix stop with \"NAs not",
  "  allowed in input\" (`waic()` too, by the same check). Seen at the",
  "  consolidation on a constructed draws object (`dev/rel081-loo-exact.R`,",
  "  `dev/rel081-loo-inf.R`). A sampled scalar `xi` does not reach it: the",
  "  sampler reads it through brms's `tmp_xi` chart, which keeps every",
  "  fitted row inside the support. Where that chart is blocked (a",
  "  predicted `xi`, a `cens()` or `trunc()` term, a second GEV response,",
  "  a mixture component) the sampler reads the fitted objective, whose",
  "  barrier is finite just outside a row's support, so a draw can reach",
  "  it there. brms rejects such a draw in its sampler."),
L("- **MINOR. `loo()`, `waic()` and `psis()` of draws stop on a `-Inf`",
  "  cell, by decision.** With the exact reader, a draw that puts a",
  "  `gen_extreme_value()` row outside its support gives `-Inf` in",
  "  `log_lik()`. At the release all three stop with one error of class",
  "  `frmtmb_loglik_not_finite` (`loo_check_finite()` in frmtmb.sample's",
  "  `loo_matrix()`), as brms 2.23.0 stops on its `NaN` there (\"NAs not",
  "  allowed in input\"). Before that, `loo()` and `psis()` stopped in the",
  "  loo package (\"All input values must be finite or -Inf\": loo 2.10.1's",
  "  `validate_ll()` refuses the `+Inf` that `-Inf` becomes in the",
  "  importance ratios) and `waic()` gave an elpd of `NaN`",
  "  (`dev/rel081-loo-inf.R`). Seen only on a draws object built by hand",
  "  (`dev/rel081-loo-exact.R`). A sampled scalar `xi` cannot reach it:",
  "  brms's `tmp_xi` chart keeps every fitted row inside the support.",
  "  Where the chart is blocked (a predicted `xi`, a `cens()` or `trunc()`",
  "  term, a second GEV response, a mixture component) the sampler reads",
  "  the fitted objective, which is continuous across the end of the",
  "  support and falls by 1e8 per unit of `t` outside it, so a draw",
  "  outside must lie in a band about 1e-8 wide: 0 of 44,000 draws in 22",
  "  such fits (release review, finding 2). Left: nothing, unless a model",
  "  is found whose draws reach the band."),
"test-backlog.md: the -Inf item")

## 5. the commit script --------------------------------------------------
cs <- "dev/rel081-commit.sh"
edit(cs,
L("# main's checkout must be clean: the script stops before any write if it",
  "# is not."),
L("# main's checkout must be clean: the script stops before any write if it",
  "# is not.",
  "#",
  "# Resuming after a partial run. The lane commits are made one at a time.",
  "# If the script stops after some of them, those lanes' branches are one",
  "# commit past e0cce906, and the HEAD checks stop a rerun. Do not reset",
  "# them. For each lane already committed, check `git -C <tree> log",
  "# --oneline -2` (its commit, then e0cce906) and `git -C <tree> status",
  "# --porcelain` (nothing but the excluded folders). Then run the",
  "# remaining commit_lane calls and the release block by hand, in this",
  "# file's order: the release commit's parents are e0cce906 and the five",
  "# wt-r81-<lane> branches, whichever run made them. If the release",
  "# commit exists and main did not move, run only",
  "# `git -C <main checkout> merge --ff-only <release commit>`."),
"commit script: how to resume")
edit(cs,
L("      { echo \"$t: no index to read\" >&2; return 1; }",
  "  done",
  "}",
  "for l in $LANES; do"),
L("      { echo \"$t: no index to read\" >&2; return 1; }",
  "  done",
  "}",
  "# read before the first git call on any tree, so that the dry run's",
  "# \"unchanged\" covers every read the script makes",
  "before=$(index_sums) || exit 1",
  "for l in $LANES; do"),
"commit script: index sums before the first git call")
edit(cs,
L("if $DRY; then",
  "  before=$(index_sums) || exit 1",
  "  mkdir -p"),
L("if $DRY; then",
  "  mkdir -p"),
"commit script: the old place of the index sums")
edit(cs,
L("         dev/rel081-release-fix2.R dev/rel081-news.R dev/rel081-backlog.R \\"),
L("         dev/rel081-release-fix2.R dev/rel081-release-fix3.R \\",
  "         dev/rel081-news.R dev/rel081-backlog.R \\"),
"commit script: fix3 among the required files")
edit(cs,
L("commit_lane() {",
  "  local n=\"$1\" subject=\"$2\" body=\"$3\"",
  "  ( cd \"$R/frmtmb-wt-$n\" && git add -A \"${SPEC[@]}\" )"),
L("# the paths now staged in a tree, as take_list spells them",
  "staged_list() {",
  "  ( cd \"$R/frmtmb-wt-$1\" &&",
  "      git diff --cached --name-only --no-renames -z | tr '\\0' '\\n' | sort -u )",
  "}",
  "# git add must stage exactly what take_list listed before it, or the",
  "# commit would hold something the dry run did not show",
  "add_checked() {",
  "  local want got",
  "  want=$(take_list \"$1\") || exit 1",
  "  ( cd \"$R/frmtmb-wt-$1\" && git add -A \"${SPEC[@]}\" )",
  "  got=$(staged_list \"$1\") || { echo \"wt-$1: git diff --cached failed\"; exit 1; }",
  "  [ \"$want\" = \"$got\" ] ||",
  "    { echo \"wt-$1: the staged paths differ from the list; nothing committed here\"; exit 1; }",
  "}",
  "commit_lane() {",
  "  local n=\"$1\" subject=\"$2\" body=\"$3\"",
  "  add_checked \"$n\""),
"commit script: lanes stage exactly take_list's paths")
edit(cs,
L("git -C \"$REL\" add -A \"${SPEC[@]}\"",
  "tree=$(git -C \"$REL\" write-tree)"),
L("add_checked rel081",
  "tree=$(git -C \"$REL\" write-tree)"),
"commit script: the release tree stages exactly take_list's paths")
edit(cs,
L("brief's merge steps), dev/rel081-release-fix2.R (the upstream bug",
  "number), dev/rel081-news.R, dev/rel081-backlog.R."),
L("brief's merge steps), dev/rel081-release-fix2.R (the upstream bug",
  "number), dev/rel081-release-fix3.R (the release review's five",
  "minors), dev/rel081-news.R, dev/rel081-backlog.R."),
"commit script: the release message names fix3")

## 5. brms entries of dev/upstream-bugs.md in numeric order --------------
ub <- file.path(REL, "dev/upstream-bugs.md")
f <- rd(ub)
x <- strsplit(f$txt, "\n", fixed = TRUE)[[1L]]
s <- which(x == "## brms")
e <- which(grepl("^## ", x) & seq_along(x) > s)[1L]
stopifnot(length(s) == 1L, !is.na(e))
h <- which(grepl("^### brms-[0-9]+\\. ", x) & seq_along(x) > s & seq_along(x) < e)
num <- as.integer(sub("^### brms-([0-9]+)\\..*", "\\1", x[h]))
if (!is.unsorted(num)) {
  cat("already applied: brms entries in numeric order\n")
} else {
  ends <- c(h[-1L] - 1L, e - 1L)
  blocks <- lapply(seq_along(h), function(j) {
    b <- x[h[j]:ends[j]]
    # each block ends with one blank line before the next heading
    while (length(b) && !nzchar(b[length(b)])) b <- b[-length(b)]
    c(b, "")
  })
  x2 <- c(x[seq_len(h[1L] - 1L)], unlist(blocks[order(num)]), x[e:length(x)])
  stopifnot(length(x2) == length(x), setequal(x2, x))
  wrt(f, paste(x2, collapse = "\n"), ub)
  cat("applied: brms entries in numeric order (", paste(num, collapse = " "),
      "->", paste(sort(num), collapse = " "), ")\n")
}

## 3 and 4: the NEWS files ----------------------------------------------
rs <- file.path(R.home("bin"), "Rscript")
st <- system2(rs, shQuote(file.path(REL, "dev/rel081-news.R")), stdout = TRUE,
              stderr = TRUE)
writeLines(st)
stopifnot(identical(st[length(st)], "NEWS DONE"))
cat("FIX3 DONE\n")
