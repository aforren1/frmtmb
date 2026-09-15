# The collision table, read from dev/samplegen-out/check-<label>/<mode>.txt
# as dev/samplegen-check.R wrote them. Counts come from the files; none
# is typed.
#   Rscript dev/samplegen-summary.R base fix
labels <- commandArgs(trailingOnly = TRUE)
modes <- c(S = "brms, then frmtmb.sample",
           T = "frmtmb.sample, then brms",
           U = "frmtmb.sample; brms only loaded",
           I = "brms; an importer loaded, probed inside",
           D = "brms, sample; detach and reattach both",
           R = "loo, sample; loo unloaded; brms loaded",
           N = "frmtmb.sample alone",
           P = "7 other owners, then frmtmb.sample",
           Q = "frmtmb.sample, then 7 other owners",
           G = "gratia, then frmtmb.sample")
field <- function(ln, key) {
  x <- grep(paste0("^", key, " "), ln, value = TRUE)
  if (!length(x)) return(NA_character_)
  trimws(sub(paste0("^", key, " "), "", x[1]))
}
brms_loaded <- c("S", "T", "U", "I", "D", "R")
cat("```\n== frmtmb.sample collision table, dev/samplegen-check.R ==\n")
cat("one process per row; lookup in the method table of environment(generic)\n")
cat("brmsfit: brms's class method unreachable (rows where brms is loaded)\n")
cat("draws:   the method that would run for a frmtmb_draws is not this\n")
cat("         package's own (a foreign .default counts as a failure)\n")
cat("owners:  owner class methods the CONTROL reaches and the test does not\n\n")
cat(sprintf("%-5s %-4s %-41s %-9s %-7s %s\n", "build", "mode", "order",
            "brmsfit", "draws", "owners"))
for (lb in labels) {
  for (m in names(modes)) {
    f <- file.path("dev/samplegen-out", paste0("check-", lb),
                   paste0(m, ".txt"))
    if (!file.exists(f)) { cat(lb, m, "MISSING\n"); next }
    ln <- readLines(f, warn = FALSE)
    if (!any(ln == "CHILDOK")) { cat(lb, m, "DID NOT FINISH\n"); next }
    b <- sub(" of ", "/", field(ln, "BRMSFIT_LOST"))
    if (!m %in% brms_loaded) b <- "n/a"
    d <- sub(" of ", "/", field(ln, "DRAWS_NOT_OURS"))
    o <- sub(" of (\\d+).*", "/\\1", field(ln, "OWNER_METHODS_LOST"))
    if (m == "I") o <- "n/a"
    cat(sprintf("%-5s %-4s %-41s %-9s %-7s %s\n", toupper(lb), m,
                modes[[m]], b, d, o))
    extra <- c(
      if (m %in% brms_loaded && nzchar(field(ln, "BRMSFIT_SILENT_DEFAULT")))
        paste("  brmsfit silently to a .default:",
              field(ln, "BRMSFIT_SILENT_DEFAULT")),
      if (nzchar(field(ln, "DRAWS_FOREIGN_DEFAULT")))
        paste("  draws to a foreign .default:",
              field(ln, "DRAWS_FOREIGN_DEFAULT")),
      if (m == "G") paste("  posterior_samples.gam reachable:",
                          field(ln, "GAM_REACHABLE")),
      if (m == "R") paste("  psis resolves to", field(ln, "PSIS_FROM"),
                          " stale:", field(ln, "PSIS_STALE")),
      paste("  resolves:", field(ln, "RESOLVES"),
            "| active in namespace:", field(ln, "ACTIVE_NS"),
            if (m == "I") paste("| in importer:",
                                field(ln, "ACTIVE_IMPORTER")) else ""))
    if (m %in% c("S", "T", "I", "G", "R", "Q")) cat(extra, sep = "\n")
  }
}
cat("```\n")
