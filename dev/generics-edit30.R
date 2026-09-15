source("C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-sub.R")

# The published ratio was 443 at the median row and 768.7 at row 1.
# Re-running against the SAME cached brms fit gives 441.3 and 765.1,
# because brms's `predict()` summarises FRESH predictive draws: it is a
# Monte Carlo quantity and four significant figures on it are three
# too many. The `fitted()` agreement figures are unchanged to the last
# digit, because both sides of that comparison are deterministic given
# the fit, and those are the ones the claim rests on.
sub1("R/scales.R",
paste0("#' by a factor of 443 at the MEDIAN row, and 768.7 at row 1, because\n",
       "#' it is a different scale rather than a different answer.\n"),
paste0("#' by roughly 440 at the median row and 765 at row 1, because it is\n",
       "#' a different scale rather than a different answer. Those two are\n",
       "#' given loosely on purpose: brms's `predict()` summarises fresh\n",
       "#' predictive draws, so it moves between calls on one fit, where\n",
       "#' the `fitted()` figures above are deterministic and repeat to the\n",
       "#' last digit.\n"))
cat("DONE\n")
