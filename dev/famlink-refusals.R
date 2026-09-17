# Every (family, roster link) pair brms 2.23.0 REFUSES for the mean, and
# whether one frmtmb library constructs it.
#
#   Rscript dev/famlink-refusals.R base
#   Rscript dev/famlink-refusals.R lane
#
# The pair list comes from R/links-brms.R (the generated brms table) and
# the link roster frmtmb registers, read from the SOURCE so both arms try
# the same list. multinomial has no link argument and is left out.
arm <- commandArgs(trailingOnly = TRUE)[1]
lib <- switch(arm,
              base = "C:/Users/adf44/source/r/rellib-r3",
              lane = "C:/Users/adf44/source/r/famlink-lib",
              stop("arm must be base or lane"))
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
tab_env <- new.env()
sys.source("R/links-brms.R", envir = tab_env)
sets <- tab_env$brms_mu_links
roster <- names(frmtmb:::frmtmb_links)
rows <- list()
for (fam in setdiff(names(sets), "multinomial")) {
  ctor <- frmtmb:::family_registry[[switch(fam, beta = "Beta", fam)]]
  for (lk in setdiff(roster, sets[[fam]])) {
    r <- tryCatch({ ctor(link = lk); "constructed" },
                  error = function(e) "refused")
    rows[[length(rows) + 1L]] <- data.frame(family = fam, link = lk,
                                            outcome = r)
  }
}
out <- do.call(rbind, rows)
cat("---- GENERATED: dev/famlink-refusals.R", arm, "----\n")
cat(sprintf("pairs brms refuses: %d over %d families\n", nrow(out),
            length(unique(out$family))))
cat(sprintf("constructed: %d; refused: %d\n", sum(out$outcome == "constructed"),
            sum(out$outcome == "refused")))
con <- out[out$outcome == "constructed", ]
if (nrow(con)) {
  byfam <- tapply(con$link, con$family, length)
  cat("constructed, by family:",
      paste0(names(byfam), " ", byfam, collapse = "; "), "\n")
}
cat("---- END GENERATED ----\n")
