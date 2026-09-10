# Lane `latent`: undo one shipped decision, on a COPY, so the test that
# pins it can be seen failing. Called by
# dev/latent-falsealarm-check.ps1.
#
#   Rscript dev/latent-falsealarm-revert.R <copy of the package> <mode>
#
# mode "notol"  round 1: the two comparisons with no tolerance at all,
#               which announced a local optimum on 6 of 6 unimodal fits.
# mode "refs"   round 2: one tolerance FUNCTION called with different
#               references at each site, which leaves a band of
#               `grad_tol` where the modes table and the printed verdict
#               disagree.
# mode "w05"    round 2: the shrink weight back at 0.5, the low edge of
#               the usable interval.
# mode "w099"   round 2: the shrink weight at 0.99, the high edge.

av <- commandArgs(trailingOnly = TRUE)
root <- av[[1L]]
mode <- if (length(av) >= 2L) av[[2L]] else "notol"
file <- if (mode %in% c("w05", "w099")) "lca.R" else "hmm-starts.R"
p <- file.path(root, "R", file)
src <- paste(readLines(p, warn = FALSE), collapse = "\n")

sub1 <- function(x, from, to) {
  if (!grepl(from, x, fixed = TRUE)) {
    stop("anchor not found: ", from, call. = FALSE)
  }
  sub(from, to, x, fixed = TRUE)
}

out <- src
if (identical(mode, "notol")) {
  out <- sub1(out,
              "    if (conv && is.finite(lli) && lli > bl + mode_tol) {",
              "    if (conv && is.finite(lli) && lli > bl) {")
  out <- sub1(out, "  if (gap > x[[\"mode_tol\"]]) {", "  if (gap > 0) {")
} else if (identical(mode, "refs")) {
  out <- sub1(out,
              "    if (conv && is.finite(lli) && lli > bl + mode_tol) {",
              paste0("    if (conv && is.finite(lli) &&\n",
                     "        lli > bl + hmm_starts_tol(grad_tol, bl)) {"))
  out <- sub1(out, "  modes <- hmm_starts_modes(conv_ll, mode_tol)",
              "  modes <- hmm_starts_modes(conv_ll, grad_tol)")
  out <- sub1(out, "hmm_starts_modes <- function(ll, mode_tol) {",
              "hmm_starts_modes <- function(ll, grad_tol) {")
  out <- sub1(out, "    same <- abs(o[i] - ref) <= mode_tol",
              "    same <- abs(o[i] - ref) <= hmm_starts_tol(grad_tol, ref)")
  out <- sub1(out, "  if (gap > x[[\"mode_tol\"]]) {",
              "  if (gap > hmm_starts_tol(x[[\"grad_tol\"]], bll)) {")
} else if (identical(mode, "w05")) {
  out <- sub1(out, "      p <- 0.1 * p + 0.9 * pool",
              "      p <- 0.5 * p + 0.5 * pool")
} else if (identical(mode, "w099")) {
  out <- sub1(out, "      p <- 0.1 * p + 0.9 * pool",
              "      p <- 0.01 * p + 0.99 * pool")
} else {
  stop("unknown mode ", mode, call. = FALSE)
}
if (identical(out, src)) stop("the revert was a no-op", call. = FALSE)
writeLines(strsplit(out, "\n", fixed = TRUE)[[1L]], p)
cat("mode", mode, "reverted in", p, "\n")
