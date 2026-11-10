# Release fixes found at the consolidation of 0.81.0, applied once after
# dev/rel081-release-fix.R:
#
#   Rscript dev/rel081-release-fix2.R
#
# 1. dev/upstream-bugs.md: lane translator filed its brms 2.23.0 GEV
#    mixture finding as brms-49, but brms-49 (predictor.bprepnl() drops
#    the ARMA term) and brms-50 were already filed. It becomes brms-51,
#    everywhere this release's records cite it (the review keeps its own
#    text).
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel081"
rd <- function(p) {
  s <- rawToChar(readBin(p, "raw", file.info(p)$size))
  list(text = gsub("\r\n", "\n", s, fixed = TRUE),
       crlf = grepl("\r\n", s, fixed = TRUE))
}
wr <- function(p, f) {
  s <- if (f$crlf) gsub("\n", "\r\n", f$text, fixed = TRUE) else f$text
  writeBin(charToRaw(s), p)
}
fix <- function(rel, old, new, n = 1L) {
  p <- file.path(REL, rel)
  f <- rd(p)
  k <- lengths(regmatches(f$text, gregexpr(old, f$text, fixed = TRUE)))
  if (k != n) stop(rel, ": the edit matched ", k, " times, not ", n, ": ",
                   substr(old, 1, 60))
  f$text <- gsub(old, new, f$text, fixed = TRUE)
  wr(p, f)
  cat("edited", rel, "(", k, ")\n")
}
fix("dev/upstream-bugs.md",
    "### brms-49. A mixture's `gen_extreme_value` component must cover every response",
    "### brms-51. A mixture's `gen_extreme_value` component must cover every response")
fix("dev/upstream-bugs.md",
    "- Record: `dev/translator-findings.md`, p021.\n- Status: not reported.",
    paste0("- Record: `dev/translator-findings.md`, p021. Filed by the lane as\n",
           "  brms-49; numbered brms-51 at the 0.81.0 release, since brms-49 and\n",
           "  brms-50 were taken.\n- Status: not reported."))
fix("dev/equiv-data/adjudication.tsv", "(brms-49, dev/upstream-bugs.md)",
    "(brms-51, dev/upstream-bugs.md)")
f <- rd(file.path(REL, "dev/translator-findings.md"))
k <- lengths(regmatches(f$text, gregexpr("brms-49", f$text, fixed = TRUE)))
fix("dev/translator-findings.md", "brms-49", "brms-51", n = k)
fix("dev/translator-findings.md",
    "# Lane translator: the survey's weak and unmeasured pairs, lifted to L1\n",
    paste0("# Lane translator: the survey's weak and unmeasured pairs, lifted to L1\n\n",
           "(At the 0.81.0 release the lane's brms-49 was numbered brms-51 in\n",
           "`dev/upstream-bugs.md`, since brms-49 and brms-50 were taken; this\n",
           "record and the adjudication row say brms-51. The review keeps its\n",
           "own text.)\n"))
cat("RELEASE FIX 2 DONE\n")
