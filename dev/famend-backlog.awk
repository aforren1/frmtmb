# Lane famend: the CLOSED lines and the filed section, inserted into
# dev/test-backlog.md by its 0.82.0 line numbers in one pass.
#   awk -f dev/famend-backlog.awk dev/test-backlog.md > new && mv new ...
BEGIN {
  P = "  CLOSED in the round of 2026-11-13 (lane famend): "
  after[2750] = P "a cs() effect that closes a row's own category is named by mixture_ord_degenerate_check(), with the row and the width; on dev/optima-csmix-wall.R's 20 data sets the 18 code-1 fits (widths 1.4e-15 to 2.3e-13) are named and the 2 converged ones (0.509, 0.893) are not."
  after[3005] = P "does not reproduce on 0.81.0 or 0.82.0: 0 of 40 fits error, from the default start or the true values (dev/famend-probe-exg.R); 0.73.0's Mills-ratio density fixed it."
  after[3066] = P "frechet_fit_check() warns (frmtmb_frechet_nu_edge) where nu is within 1e-3 of 1 and the log-likelihood does not fall along the scale-held ridge toward it; seed 9031 warns, and 77 interior maxima of the 80 probe fits do not (dev/famend-probe-frechet.R)."
  after[3074] = P "the GEV's cap on exp(-w) is exp(500), which bounds the derivatives; 0 of 6000 grid points have a non-finite gradient at xi 0.06, 0.3 and 1, against 80 of 6000 (dev/famend-repro.R, gevgrad)."
  after[3081] = P "gev_fit_check() gives the message frmtmb_gev_xi_nonregular below -0.5, with confint(method = \"profile\"); fits at xi -0.663, -0.736 and -0.936 get it, -0.313 and -0.473 do not."
  after[3563] = P "a start where a prior's density is 0 moves inside its support (prior_start_inside()); the fit reaches the objective of the lb = 0 fit, which also stuck before (dev/famend-probe-gprior.R)."
  after[3619] = P "residuals(type = \"osa\") reads a censored count on the window (lo, hi]; poisson, negbinomial, binomial and zero-inflated poisson match the exact randomized PIT to 1.5e-13 (dev/famend-osa.R)."
  after[3648] = P "fixed by lane boundfam at 0.82.0: both spellings name the limit, and their log-likelihoods agree to 1e-8 relative on 20 of 20 ulp-perturbed data sets (dev/famend-perturb.R)."
  after[5725] = P "measured and not adopted: RTMB::pbeta() is 12.6 times cheaper per sweep, but on the one-step tape of a binomial (1 | g) fit it gave 704 of 800 non-finite residuals, and its log tail is -Inf at size 2000 (dev/famend-cost.R)."
  after[5728] = P "a continuous response read off its CDF gives the fresh-copy residual for a weight in (0, 1); gaussian rows match qnorm(F(y)) to 2.3e-15, and a count keeps NA (dev/famend-osa.R)."
  before[4796] = "## Filed in the round of 2026-11-13\n\n### From lane famend\n\n" \
    "- LOW. **A right-censored poisson under single-block quadrature still\n" \
    "  costs about 2 MB a row** (dev/famend-quadmem.R, the release review's\n" \
    "  construction at 1440 rows): 2773 MB and 33.4 s after the Horner upper\n" \
    "  sum of this round, against 3380 MB and 64.3 s before and 495 MB at\n" \
    "  0.72.0. The lower sum and the blend of `poisson_cdf_log()` remain,\n" \
    "  and the gamma, negative binomial and incomplete beta pairs keep\n" \
    "  their 140 to 226 times `ppois()`. The residue of the CDF-pair cost\n" \
    "  item of lane cdfs.\n\n"
}
{
  if (NR in before) printf "%s", before[NR]
  print
  if (NR in after) print after[NR]
}
