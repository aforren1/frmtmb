LIB <- "C:/Users/adf44/source/r/wt-records-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-records"

out <- file.path(WT, "dev", "records-me-rd.txt")
tools::Rd2txt(file.path(WT, "man", "frmtmb-me.Rd"), out = out)
txt <- readLines(out, warn = FALSE)

# The Value section is what pkgcheck reported missing, so grep the
# RENDERED text rather than the Rd source.
i <- grep("^Value:", txt)
cat("Value headings found:", length(i), "\n")
if (length(i)) {
  cat(paste(txt[i[1]:min(length(txt), i[1] + 12)], collapse = "\n"), "\n")
}

# the ad-env.R fix must leave no Rd behind for the internal object
cat("nl_rtmb_shadow Rd present:",
    file.exists(file.path(WT, "man", "nl_rtmb_shadow.Rd")), "\n")

# every man page still parses
rds <- list.files(file.path(WT, "man"), pattern = "[.]Rd$", full.names = TRUE)
bad <- character(0)
for (f in rds) {
  ok <- tryCatch({ tools::parse_Rd(f); TRUE },
                 error = function(e) { bad <<- c(bad, basename(f)); FALSE })
}
cat("Rd files parsed:", length(rds), " failures:", length(bad), "\n")
if (length(bad)) cat(paste(bad, collapse = ", "), "\n")

# and no Rd lacks \value that has \usage (the pkgcheck rule's shape)
novalue <- character(0)
for (f in rds) {
  p <- tools::parse_Rd(f)
  tg <- vapply(p, function(x) attr(x, "Rd_tag"), character(1))
  if ("\\usage" %in% tg || "\\name" %in% tg) {
    if (!("\\value" %in% tg) && !("\\docType" %in% tg)) {
      novalue <- c(novalue, basename(f))
    }
  }
}
cat("Rd with no \\value:", length(novalue), "\n")
if (length(novalue)) cat(paste(novalue, collapse = "\n"), "\n")
cat("DONE\n")
