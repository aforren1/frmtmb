# Second dump of brms 2.23.0 internals for lane ceplot: the new-level
# machinery behind sample_new_levels, and the hypothesis plot helpers.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
stopifnot(packageVersion("brms") == "2.23.0")
ns <- asNamespace("brms")
out <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-brmssrc"
nms <- ls(ns, all.names = TRUE)
sel <- grep("new_r|old_lev|new_lev|rsamples|rdraws|^expand_|nsamples|warmup",
            nms, value = TRUE)
for (nm in sel) {
  f <- get(nm, envir = ns)
  if (!is.function(f)) next
  writeLines(c(paste0(nm, " <- "), deparse(f, control = "useSource")),
             file.path(out, paste0(gsub("[^A-Za-z0-9_.]", "_", nm), ".R")))
}
cat(sort(sel), sep = "\n")
