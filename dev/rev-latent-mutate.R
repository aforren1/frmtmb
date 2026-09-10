# Reviewer: apply ONE mutation to a COPY of frmtmb.latent's
# hmm-starts.R, and refuse loudly if the text it expects is not there.
# A mutation that silently did nothing would make the suite look
# stronger than it is, which is the failure this whole exercise is
# about.
#
#   Rscript dev/rev-latent-mutate.R <copy dir> <M1|M2|M3|M4>

args <- commandArgs(trailingOnly = TRUE)
dir <- args[[1L]]
m <- args[[2L]]
p <- file.path(dir, "R", "hmm-starts.R")
txt <- readLines(p, warn = FALSE)
src <- paste(txt, collapse = "\n")

NL <- "
"

sub1 <- function(s, from, to) {
  n <- length(gregexpr(from, s, fixed = TRUE)[[1L]])
  if (!nzchar(from) || regexpr(from, s, fixed = TRUE) < 0) {
    stop("mutation target not found: ", from, call. = FALSE)
  }
  if (n != 1L) stop("mutation target is not unique: ", from, call. = FALSE)
  sub(from, to, s, fixed = TRUE)
}

out <- switch(
  m,
  ## the jitter never reaches the start: every refit begins at the
  ## incumbent's own estimates
  M1 = sub1(src,
            "    st[[cp]][] <- v[idx]",
            "    st[[cp]][] <- st[[cp]][]"),
  ## a refit that did not converge may be returned as $best
  M2 = sub1(src,
            "    if (conv && is.finite(lli) && lli > as.numeric(stats::logLik(best))) {",
            "    if (is.finite(lli) && lli > as.numeric(stats::logLik(best))) {"),
  ## the converged spread is taken over everything that finished
  M3 = sub1(src,
            '  conv_ll <- rows$logLik[rows$status %in% c("original", "converged")]',
            '  conv_ll <- rows$logLik[rows$status != "error"]'),
  ## the incumbent is exempt from the test the refits are held to
  M4 = sub1(src,
            "  original_converged <- is.finite(rows$grad_rel[1L]) &&\n    rows$grad_rel[1L] <= grad_tol",
            "  original_converged <- TRUE"),
  ## FIX candidate, not a defect: give the two untoleranced
  ## comparisons the same relative tolerance hmm_starts_modes()
  ## already uses, so the summary stops calling optimizer noise a
  ## local optimum
  F1 = sub1(
    sub1(src,
         "    if (conv && is.finite(lli) && lli > as.numeric(stats::logLik(best))) {",
         paste0("    .bl <- as.numeric(stats::logLik(best))", NL,
                "    if (conv && is.finite(lli) &&", NL,
                "        lli > .bl + grad_tol^2 * max(abs(.bl), 1)) {")),
    "  if (gap > 0) {",
    "  if (gap > x[[\"grad_tol\"]]^2 * max(abs(bll), 1)) {"),
  stop("unknown mutant ", m, call. = FALSE))

if (identical(out, src)) stop("mutation was a no-op", call. = FALSE)
writeLines(strsplit(out, "\n", fixed = TRUE)[[1L]], p)
cat("mutant", m, "applied to", p, "\n")
