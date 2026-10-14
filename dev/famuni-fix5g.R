p <- "extensions/frmtmb.sample/R/sample.R"
x <- paste(readLines(p), collapse = "\n")
r <- function(old, new) {
  if (!grepl(old, x, fixed = TRUE)) stop("not found: ", substr(old, 1, 60))
  x <<- sub(old, new, x, fixed = TRUE)
}
r("                          if (isTRUE(s$natural)) \"  [natural scale]\"",
  "                          # a scalar xi's scale has a note of its own
                          if (identical(s$class, \"xi\")) \"\"
                          else if (isTRUE(s$natural)) \"  [natural scale]\"")
r("    for (tm in gev_xi_chart_terms(fit)) {
      if (!identical(tm$resp, rspec$resp_name) || !length(tm$blocked)) next
      notes <- c(notes, gev_xi_blocked_text(tm, of))
    }",
"    for (tm in gev_xi_chart_terms(fit)) {
      if (!identical(tm$resp, rspec$resp_name)) next
      notes <- c(notes, if (length(tm$blocked)) {
        gev_xi_blocked_text(tm, of)
      } else {
        paste0(\"xi\", of, \" is sampled as brms samples a scalar xi, \",
               \"through tmp_xi, mapped into the support each draw's mu and \",
               \"sigma leave, and above -1, where brms's range goes on; a \",
               \"prior on class xi is a prior on tmp_xi, as in brms\")
      })
    }")
writeLines(strsplit(x, "\n")[[1]], p)
