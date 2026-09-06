# brms's count of families per link, against what core's registry has
# and what core's family roster can actually reach.
ns <- asNamespace("brms")
fns <- setdiff(grep("^[.]family_", ls(ns, all.names = TRUE), value = TRUE),
               ".family_custom")
famlinks <- list()
for (f in fns) {
  info <- tryCatch(get(f, envir = ns)(), error = function(e) NULL)
  if (is.null(info) || !length(info$links)) next
  famlinks[[sub("^[.]family_", "", f)]] <- info$links
}
allk <- sort(unique(unlist(famlinks)))
core <- names(frmtmb:::frmtmb_links)
reg <- frmtmb:::family_registry
cat(sprintf("%-14s %8s %8s %10s %s\n", "link", "brms fam",
            "in core", "core fam", "brms families core lacks"))
for (k in allk) {
  fams <- names(famlinks)[vapply(famlinks, function(v) k %in% v, TRUE)]
  alias <- c(gamma = "Gamma", beta = "Beta")
  cn <- function(f) if (f %in% names(alias)) alias[[f]] else f
  have <- vapply(fams, function(f) !is.null(reg[[cn(f)]]), TRUE)
  cat(sprintf("%-14s %8d %8s %10d %s\n", k, length(fams),
              if (k %in% core) "yes" else "NO", sum(have),
              paste(fams[!have], collapse = " ")))
}
cat("\nregistry entries brms has no name for: ",
    paste(setdiff(core, allk), collapse = ", "), "\n", sep = "")
cat("brms link names still absent from core: ",
    paste(setdiff(allk, core), collapse = ", "), "\n", sep = "")
