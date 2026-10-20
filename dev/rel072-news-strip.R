# The user's decision after the release review (option b): the 0.72.0
# NEWS sections name no lane and no `dev/` path. Those mean nothing to
# a user, and `dev/` is not in the built package. The trail stays in
# dev/round-20261020.md ("Where each NEWS bullet came from", from
# dev/rel072-news-map.md, which dev/rel072-news-merge.R writes).
#
# dev/rel072-news-merge.R sources this after the moves and the rewrap.
# Each edit is a fixed text, matched once in the file's top section
# with any run of whitespace standing for a space, so it does not
# depend on where the rewrap broke the lines; the section is rewrapped
# again after. The consolidation's own leads (dev/rel072-news-*-lead.md)
# name no lane either. `mentions()` counts what is left: a lane name,
# the word "lane", or a `dev/` path.
strip_lanes <- c("testiso", "catmulti", "prilabel", "sampleprior",
                 "bin2port", "ndaterm", "illse2", "adhess", "drop1raw")
mentions <- function(x) {
  pat <- paste0("\\blanes?\\b|dev/|\\b(",
                paste(strip_lanes, collapse = "|"), ")\\b")
  sum(vapply(gregexpr(pat, x, ignore.case = TRUE, perl = TRUE),
             function(m) sum(m > 0), 0))
}
strip_edits <- list(
  "NEWS.md" = list(
    c(" (`dev/drop1raw-findings.md`). brms has no `drop1()`.",
      ". brms has no `drop1()`."),
    c(paste(" No package code changed (lane `bin2port`,",
            "`dev/bin2port-findings.md`)."),
      " No package code changed.")),
  "extensions/frmtmb.sample/NEWS.md" = list(
    c(paste(" The table of every brms family's defaults against",
            "frmtmb's is in `dev/sampleprior-findings.md`."),
      ""))
)
# the edit, as a regex over the top section: whitespace runs match any
# whitespace, every other character literally
ws_regex <- function(s) {
  parts <- strsplit(s, "\\s+")[[1]]
  esc <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", parts)
  paste(esc, collapse = "\\s+")
}
