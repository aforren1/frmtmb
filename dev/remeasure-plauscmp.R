# The estimate comparison with brms, build against build: per model, the
# frmtmb optimizer code and the median and max abs z, so a fit that
# stopped converging or moved is named, not hidden in a median.
#
#   Rscript dev/remeasure-plauscmp.R <tag> [<tag> ...]
#
# Reads dev/remeasure-out/<tag>/port/plausibility.rds (written by
# dev/brms-port/plausibility.R through dev/remeasure-run.sh). Every tag
# is compared against the same brms reference fits, so a difference is
# frmtmb's.
root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
tags <- commandArgs(trailingOnly = TRUE)
P <- lapply(tags, function(t) {
  readRDS(file.path(root, "remeasure-out", t, "port", "plausibility.rds"))
})
names(P) <- tags
ids <- P[[1]]$table$id
f3 <- function(x) ifelse(is.na(x), "-", formatC(x, digits = 3, format = "fg"))
col <- function(t, v) P[[t]]$table[[v]][match(ids, P[[t]]$table$id)]
cat("## per model: frmtmb status / optimizer code / median abs z / max",
    "abs z\n\n")
cat(sprintf("%-26s %s\n", "model", paste(sprintf("%-24s", tags),
                                          collapse = " ")))
changed <- character()
for (k in seq_along(ids)) {
  cells <- vapply(tags, function(t) {
    sprintf("%-5s %2s %6s %6s", col(t, "frmtmb")[k],
            ifelse(is.na(col(t, "fconv")[k]), "-",
                   as.character(col(t, "fconv")[k])),
            f3(col(t, "med_z")[k]), f3(col(t, "max_z")[k]))
  }, "")
  cat(sprintf("%-26s %s\n", sub("^brms_", "", ids[k]),
              paste(sprintf("%-24s", cells), collapse = " ")))
  st <- vapply(tags, function(t) col(t, "frmtmb")[k], "")
  cv <- vapply(tags, function(t) as.character(col(t, "fconv")[k]), "")
  mz <- vapply(tags, function(t) col(t, "max_z")[k], 0)
  if (length(unique(st)) > 1 || length(unique(cv)) > 1 ||
      (all(!is.na(mz)) && diff(range(mz)) > 0.01)) {
    changed <- c(changed, ids[k])
  }
}
cat("\nmodels whose status, optimizer code or max abs z (by more than",
    "0.01) differ between builds:", length(changed), "\n")
for (id in changed) {
  cat("\n###", id, "\n")
  for (t in tags) {
    d <- P[[t]]$detail[[id]]
    w <- P[[t]]$table$fwarn[P[[t]]$table$id == id]
    cat(sprintf("  %-6s code %s %s\n", t,
                P[[t]]$table$fconv[P[[t]]$table$id == id], w))
    if (!is.null(d$fx)) {
      cat("        ", paste(sprintf("%s=%s(z %s)", d$fx$name, f3(d$fx$frmtmb),
                                     f3(d$fx$z)), collapse = " "), "\n")
    }
  }
}
both_ok <- function(t) {
  R <- P[[t]]$table
  R$id[R$frmtmb == "OK" & R$brms == "OK" & R$n_b > 0 & !is.na(R$rhat) &
         R$rhat <= 1.05 & !is.na(R$fconv) & R$fconv == 0]
}
zs <- function(t, ids) unlist(lapply(P[[t]]$detail[ids],
                                     function(d) abs(d$fx$z)))
line <- function(t, ids) {
  z <- zs(t, ids)
  cat(sprintf(paste("%-6s %2d models, %3d coefficients; abs z median %s,",
                    "90th pct %s, max %s; > 1: %d\n"), t, length(ids),
              length(z), f3(stats::median(z)), f3(stats::quantile(z, 0.9)),
              f3(max(z)), sum(z > 1)))
}
cat("\n## both converged, per build (each build's own set)\n\n")
for (t in tags) line(t, both_ok(t))
# Two builds' own sets differ when a fit starts or stops converging, so
# their medians are not comparable; the common set is.
common <- Reduce(intersect, lapply(tags, both_ok))
cat("\n## both converged on EVERY build: the common set\n\n")
for (t in tags) {
  cat(sprintf("%-6s not in the common set: %s\n", t,
              paste(setdiff(both_ok(t), common), collapse = ", ")))
}
for (t in tags) line(t, common)
cat("\n## end\n")
