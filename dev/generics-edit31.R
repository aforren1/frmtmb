source("C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-sub.R")

# A REGRESSION found by `R CMD check`, which runs the whole suite in
# ONE process, and which the lane's own one-file-per-process run could
# not see. `hypothesis()`'s GENERIC did work: it armed the
# reserved-name shadowing note around `UseMethod()`. As soon as
# frmtmb's binding resolves to brms's generic, that work is gone and
# the note never fires. test-naming-collisions.R asserts the note and
# failed 8 assertions, but only after some earlier file had loaded
# brms's namespace, which is why one file per process passed.
#
# Both round-1 designs had this; it is not new to the active binding.
# The fix is to stop putting logic in a generic that is meant to be
# replaceable: arm and disarm inside the METHODS, where it runs
# whichever generic dispatched. The save-and-restore in
# hyp_shadow_arm()/disarm() already makes nesting safe.

sub1("R/confint.R",
paste0("hypothesis <- function(x, ...) {\n",
       "  # the shadowing note belongs to the call the user typed, not to any of\n",
       "  # the many environment rebuilds it triggers, so it is armed here and\n",
       "  # restored when the method returns (on.exit survives UseMethod)\n",
       "  old <- hyp_shadow_arm()\n",
       "  on.exit(hyp_shadow_disarm(old), add = TRUE)\n",
       "  UseMethod(\"hypothesis\")\n",
       "}\n"),
paste0("hypothesis <- function(x, ...) UseMethod(\"hypothesis\")\n"))

sub1("R/confint.R",
paste0("hypothesis.frmtmb_fit <- function(x, hypothesis, alpha = 0.05,\n",
       "                                  method = c(\"wald\", \"profile\", \"boot\"),\n",
       "                                  nsim = 500, seed = NULL, class = NULL,\n",
       "                                  group = NULL, vcov = NULL, ...) {\n",
       "  method <- match.arg(method)\n"),
paste0("hypothesis.frmtmb_fit <- function(x, hypothesis, alpha = 0.05,\n",
       "                                  method = c(\"wald\", \"profile\", \"boot\"),\n",
       "                                  nsim = 500, seed = NULL, class = NULL,\n",
       "                                  group = NULL, vcov = NULL, ...) {\n",
       "  # The shadowing note belongs to the call the user typed, not to\n",
       "  # any of the many environment rebuilds it triggers, so it is\n",
       "  # armed once here and restored when this method returns.\n",
       "  #\n",
       "  # It is armed in the METHOD and not in the generic on purpose.\n",
       "  # frmtmb's exported `hypothesis` is an active binding that\n",
       "  # resolves to brms's generic whenever brms is loaded\n",
       "  # (R/generic-owners.R), and brms's generic is a bare\n",
       "  # UseMethod(): anything this package puts in its own generic is\n",
       "  # simply not run then. Measured, when it was in the generic:\n",
       "  # test-naming-collisions.R lost 8 assertions under R CMD check,\n",
       "  # which runs the suite in one process where an earlier file had\n",
       "  # already loaded brms.\n",
       "  old <- hyp_shadow_arm()\n",
       "  on.exit(hyp_shadow_disarm(old), add = TRUE)\n",
       "  method <- match.arg(method)\n"))

sub1("R/multiple.R",
paste0("hypothesis.frmtmb_multiple <- function(x, hypothesis, alpha = 0.05,\n",
       "                                       class = NULL, group = NULL, ...) {\n",
       "  if (...length()) {\n"),
paste0("hypothesis.frmtmb_multiple <- function(x, hypothesis, alpha = 0.05,\n",
       "                                       class = NULL, group = NULL, ...) {\n",
       "  # armed in the method, not the generic: see hypothesis.frmtmb_fit\n",
       "  old <- hyp_shadow_arm()\n",
       "  on.exit(hyp_shadow_disarm(old), add = TRUE)\n",
       "  if (...length()) {\n"))
cat("DONE\n")
