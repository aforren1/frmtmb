setwd("C:/Users/adf44/source/r/frmtmb-wt-famuni")
ed <- function(p, old, new) {
  x <- paste(readLines(p), collapse = "\n")
  if (!grepl(old, x, fixed = TRUE)) stop("not found in ", p)
  writeLines(strsplit(sub(old, new, x, fixed = TRUE), "\n")[[1]], p)
}
ed("NEWS.md",
"    holds `xi` above -1. A fit that stops at the edge, `xi` within
    1e-3 of -1, warns that its log-likelihood is a supremum.",
"    holds `xi` above -1. A fit that stops within 1e-3 of -1 with the
    likelihood still rising toward the edge warns that its
    log-likelihood is a supremum. A maximum inside that band gets a
    message instead.")
ed("R/families-extreme.R",
"#' A fit whose `xi` stops within 1e-3 of -1 warns: the likelihood rises
#' toward the edge, so its log-likelihood is a supremum, not a maximum.
#' For `-1 < xi < -0.5` the maximum exists but is not regular (Smith,
#' 1985), so its standard errors do not have their usual meaning; no
#' message says so.",
"#' A fit whose `xi` stops within 1e-3 of -1 with the likelihood still
#' rising toward the edge warns: its log-likelihood is a supremum, not a
#' maximum. A maximum inside that band gets a message. For
#' `-1 < xi < -0.5` the maximum exists but is not regular (Smith, 1985),
#' so its standard errors do not have their usual meaning; outside the
#' band no message says so.")
cat("ok\n")
