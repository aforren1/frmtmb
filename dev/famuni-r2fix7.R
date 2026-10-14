p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
r <- function(old, new) {
  if (!grepl(old, x, fixed = TRUE)) stop("not found: ", substr(old, 1, 60))
  x <<- sub(old, new, x, fixed = TRUE)
}
r("  } else if (abs(dll) > 0.1 || dll < 0) {", "  } else if (abs(dll) > 0.1) {")
r("                \"edge, not a maximum, stopped the fit. Below -1 \",
                \"the likelihood has no maximum, so the fit stopped on \",
                \"that edge, where the log-likelihood is a supremum and \",
                \"the standard errors do not have their usual meaning. A \",",
  "                \"edge, not a maximum, stopped the fit. Below -1 \",
                \"the likelihood has no maximum; on the edge the \",
                \"log-likelihood is a supremum and the standard errors do \",
                \"not have their usual meaning. A \",")
writeLines(strsplit(x, "\n")[[1]], p)
