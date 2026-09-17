# Reviewer, lane wt-conditions (recheck 1): run the Rd examples of all 8
# packages from the INSTALLED help of one library arm, one Rscript per
# topic, and record whether each topic ran to the end or which error it
# stopped at. Base and lane are then compared by topic.
#   Rscript dev/conditions-rev-rdex.R base|lane <P>
av <- commandArgs(TRUE); arm <- av[1L]; P <- as.integer(av[2L])
root <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
lib <- if (arm == "lane") "C:/Users/adf44/source/r/conditions-lib" else
  "C:/Users/adf44/source/r/rellib-r3"
out <- file.path(root, "dev/conditions-rev-log/rdex", arm)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rscript <- "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
pkgs <- c("frmtmb", "frmtmb.coupling", "frmtmb.eam", "frmtmb.latent",
          "frmtmb.learn", "frmtmb.ode", "frmtmb.sample", "frmtmb.spline")
jobs <- list()
for (p in pkgs) {
  db <- tools::Rd_db(p, lib.loc = lib)
  for (nm in names(db)) {
    ex <- file.path(out, paste0(p, "__", sub("[.]Rd$", "", nm), ".R"))
    tools::Rd2ex(db[[nm]], ex, commentDontrun = TRUE,
                 commentDonttest = FALSE)
    if (!file.exists(ex)) next
    body <- readLines(ex, warn = FALSE)
    writeLines(c(
      sprintf(".libPaths(c('%s', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))", lib),
      "Sys.setenv(NOT_CRAN = 'true')",
      sprintf("suppressPackageStartupMessages(library(%s))", p),
      "options(warn = 1)", "set.seed(1)",
      "res <- tryCatch({", body, "'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))",
      "cat('\\nRDEX-RESULT', gsub('\\n', ' ', res), '\\n')"), ex)
    jobs[[length(jobs) + 1L]] <- ex
  }
}
cat("topics with examples:", length(jobs), "\n")
run <- function(ex) {
  o <- suppressWarnings(system2(rscript, shQuote(ex), stdout = TRUE,
                                stderr = TRUE))
  writeLines(o, sub("[.]R$", ".log", ex))
  r <- grep("^RDEX-RESULT", o, value = TRUE)
  if (!length(r)) r <- paste("RDEX-RESULT NO RESULT", tail(o, 1))
  paste(basename(ex), r)
}
cl <- parallel::makeCluster(P)
parallel::clusterExport(cl, "rscript")
res <- parallel::parLapplyLB(cl, jobs, run)
parallel::stopCluster(cl)
writeLines(sort(unlist(res)), file.path(out, "results.txt"))
cat("done\n")
