# The short check of the release fixes (dev/reviews/2026-11-10-release.md,
# "Short check of the release fixes", verdict FIX on one minor), applied
# to the release tree. Rerunnable, as dev/rel081-release-fix3.R: each edit
# is applied when its old text is there once, skipped when its new text
# is already there, and stops the script otherwise.
#
#   Rscript dev/rel081-release-fix4.R
#
# 1. pp_check()'s loo_* types build their log_lik() matrix in
#    draws_ppc_psis(), not loo_matrix(), so they stopped in the loo package
#    unclassed where NEWS and ?sample-loo say they stop like loo(). They
#    now call loo_check_finite() too, and the new test reads all three.
# 2. No BREAKING mark on the compat bullet: frm() already refused each of
#    those models, so only frm_compat() rows change (dev/rel081-news.R,
#    which this script runs).
# The commit script names this file among the required ones and in the
# release message.
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
  # an insertion keeps its old text inside the new, so it is done once
  # the new text is there; a deletion's new text is inside the old, so it
  # is done once the old text is gone
  inserted <- grepl(old, new, fixed = TRUE)
  done <- if (inserted) n_of(f$txt, new) >= 1L else n_of(f$txt, old) == 0L &&
    n_of(f$txt, new) >= 1L
  if (done) {
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

## 1. pp_check()'s loo_* types ------------------------------------------
edit("extensions/frmtmb.sample/R/methods-draws.R",
L("             \"and log_lik() refused: \", conditionMessage(e), call. = FALSE)",
  "  })",
  "  list(ll = ll, chain_id = attr(ll, \"chain_id\"))"),
L("             \"and log_lik() refused: \", conditionMessage(e), call. = FALSE)",
  "  })",
  "  # this matrix does not come from loo_matrix(), so it is checked here,",
  "  # or a cell that is not finite would stop in the loo package unclassed",
  "  loo_check_finite(ll, paste0(\"pp_check(type = '\", type, \"')\"))",
  "  list(ll = ll, chain_id = attr(ll, \"chain_id\"))"),
"methods-draws.R: draws_ppc_psis() calls loo_check_finite()")
edit("extensions/frmtmb.sample/tests/testthat/test-loo-newdata-fitted.R",
L("  expect_s3_class(suppressWarnings(loo(cs$ds)), \"loo\")",
  "})"),
L("  expect_s3_class(suppressWarnings(loo(cs$ds)), \"loo\")",
  "  # pp_check()'s loo_* types build the matrix on their own path",
  "  skip_if_not_installed(\"bayesplot\")",
  "  for (ty in c(\"loo_pit_overlay\", \"loo_intervals\", \"loo_pit_qq\")) {",
  "    expect_error(suppressWarnings(pp_check(cs$moved, type = ty)),",
  "                 class = \"frmtmb_loglik_not_finite\", info = ty)",
  "  }",
  "})"),
"test-loo-newdata-fitted.R: the three loo_* types of pp_check()")

## 2. no BREAKING mark on the compat bullet -----------------------------
edit("dev/rel081-news.R",
L("## BREAKING marks the release review asked for (finding 3): code8's",
  "## bullets 2 to 4 change standard errors and the gradient warning, and",
  "## compat's bullet changes rows of frm_compat() from \"works\" to",
  "## \"refused\". The mark lengthens the first line, so the bullet is",
  "## wrapped again at 72 columns."),
L("## BREAKING marks the release review asked for (finding 3): code8's",
  "## bullets 2 to 4 change standard errors and the gradient warning. The",
  "## mark lengthens the first line, so the bullet is wrapped again at 72",
  "## columns. compat's bullet has none: frm() already refused each of",
  "## those models, and only rows of frm_compat() change (the review's",
  "## short check)."),
"rel081-news.R: the comment on the marks")
edit("dev/rel081-news.R",
L("stopifnot(length(comp) == 1L)",
  "comp[[1L]] <- breaking(comp[[1L]])"),
L("stopifnot(length(comp) == 1L)"),
"rel081-news.R: the compat bullet keeps its own first line")

## the commit script ----------------------------------------------------
cs <- "dev/rel081-commit.sh"
edit(cs,
L("         dev/rel081-release-fix2.R dev/rel081-release-fix3.R \\"),
L("         dev/rel081-release-fix2.R dev/rel081-release-fix3.R \\",
  "         dev/rel081-release-fix4.R \\"),
"commit script: fix4 among the required files")
edit(cs,
L("number), dev/rel081-release-fix3.R (the release review's five",
  "minors), dev/rel081-news.R, dev/rel081-backlog.R."),
L("number), dev/rel081-release-fix3.R (the release review's five",
  "minors), dev/rel081-release-fix4.R (its short check),",
  "dev/rel081-news.R, dev/rel081-backlog.R."),
"commit script: the release message names fix4")

## the NEWS files -------------------------------------------------------
rs <- file.path(R.home("bin"), "Rscript")
st <- system2(rs, shQuote(file.path(REL, "dev/rel081-news.R")), stdout = TRUE,
              stderr = TRUE)
writeLines(st)
stopifnot(identical(st[length(st)], "NEWS DONE"))
cat("FIX4 DONE\n")
