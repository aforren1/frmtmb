# Reviewer: plant one hazard read in a COPY of hmm-starts.R, and refuse
# loudly if the anchor text is not there. `est$b` on a container named
# `est` partial-matches `beta` on a parameter list that has no `b`,
# which is the failure frmtmb's rule exists for.
#
#   Rscript dev/rev-latent-plant.R <copy dir>

dir <- commandArgs(trailingOnly = TRUE)[[1L]]
p <- file.path(dir, "R", "hmm-starts.R")
txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
anchor <- "hmm_starts_scale <- function(fit) {\n  par <- "
if (!grepl(anchor, txt, fixed = TRUE)) {
  stop("anchor not found; the plant would have been a no-op", call. = FALSE)
}
planted <- paste0("hmm_starts_scale <- function(fit) {\n",
                  "  est <- fit[[\"estimates\"]]\n",
                  "  if (FALSE) est$b\n  par <- ")
txt <- sub(anchor, planted, txt, fixed = TRUE)
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1L]], p)
cat("hazard planted in", p, "\n")
