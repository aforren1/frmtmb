# Lane `latent`: prove the hazard guard REACHES the new R/ file.
#
# test-bracket-access.R asserts `frm_hazard_reads("frmtmb.latent")` is
# empty. An empty result is only evidence if the scanner would have
# found a hit in the new file, so this script plants one there, and
# takes it out again. "No test reaches it" is not evidence of
# unreachability, and neither is a guard that passed.
#
# The planted read is `est$b` on a container named `est`, which is one
# of frmtmb's hazard containers: `$b` partial-matches `beta` on a
# parameter list that has no `b`, which is the exact failure the rule
# exists for.
#
#   Rscript dev/latent-hazard-check.R absent   # plant the hazard
#   Rscript dev/latent-hazard-check.R present  # take it out

source("dev/latent-env.R")
mode <- commandArgs(trailingOnly = TRUE)[[1L]]
p <- "extensions/frmtmb.latent/R/hmm-starts.R"
txt <- paste(readLines(p, warn = FALSE), collapse = "\n")

clean <- "hmm_starts_scale <- function(fit) {\n  par <- "
planted <- paste0("hmm_starts_scale <- function(fit) {\n",
                  "  est <- fit[[\"estimates\"]]\n",
                  "  if (FALSE) est$b\n  par <- ")
has_clean <- grepl(clean, txt, fixed = TRUE)
has_planted <- grepl(planted, txt, fixed = TRUE)
stopifnot(xor(has_clean, has_planted))

txt <- if (identical(mode, "absent")) {
  sub(clean, planted, txt, fixed = TRUE)
} else {
  sub(planted, clean, txt, fixed = TRUE)
}
writeLines(txt, p)
cat("mode ", mode, ": hazard read planted = ",
    grepl(planted, txt, fixed = TRUE), "\n", sep = "")
