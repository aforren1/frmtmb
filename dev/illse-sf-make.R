# Lane illse: make the SE-check batteries of the 0.70.0 consolidation
# (dev/rel070-seflags-*.R, built from the seflags review's scripts) run
# here: copies named dev/illse-sf-*.R whose "lane" arm is the lane
# library over rellib-r8 and whose "base" arm is rellib-r8 (0.70.0).
#   Rscript dev/illse-sf-make.R
d <- "C:/Users/adf44/source/r/frmtmb-wt-illse/dev"
src <- list.files(d, "^rel070-seflags-.*[.]R$", full.names = TRUE)
src <- src[!grepl("bitunits|ablate|attrib|capwarn|falarm-cmp|fx1|interplay|shortscan",
                  src)]
for (f in src) {
  s <- readLines(f)
  s <- gsub("C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel070-seflags-",
            "C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-sf-", s,
            fixed = TRUE)
  if (grepl("common[.]R$", f)) {
    s <- sub('  lane = c("C:/Users/adf44/source/r/rellib-r8",',
             '  lane = c("C:/Users/adf44/source/r/wt-illse-lib",', s,
             fixed = TRUE)
    s <- sub('           "C:/Users/adf44/source/r/rellib-r7"),',
             '           "C:/Users/adf44/source/r/rellib-r8"),', s,
             fixed = TRUE)
    s <- sub('  base = "C:/Users/adf44/source/r/rellib-r7",',
             '  base = "C:/Users/adf44/source/r/rellib-r8",', s, fixed = TRUE)
  }
  if (grepl("b1ref", f)) {
    s <- sub('.libPaths(c("C:/Users/adf44/source/r/rellib-r8",',
             paste0('.libPaths(c(if (Sys.getenv("ILLSE_ARM", "lane") == ',
                    '"lane") "C:/Users/adf44/source/r/wt-illse-lib",'),
             s, fixed = TRUE)
    s <- sub('            "C:/Users/adf44/source/r/rellib-r7",',
             '            "C:/Users/adf44/source/r/rellib-r8",', s,
             fixed = TRUE)
  }
  s <- c(paste0("# Lane illse: copy of dev/", basename(f),
                " by dev/illse-sf-make.R (lane arm: the lane library over"),
         "# rellib-r8; base arm: rellib-r8).", s)
  writeLines(s, file.path(d, sub("^rel070-seflags-", "illse-sf-",
                                 basename(f))))
}
cat(length(src), "copies\n")
