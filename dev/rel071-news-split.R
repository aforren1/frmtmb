# The release review (section 10): NEWS sentences over 25 words, split
# after the NEWS merge (dev/rel071-news-merge.R sources this), so the
# lanes' NEWS stays as each lane wrote it. Each edit is a fixed string
# that must match once in the file it names. dev/rel071-news-long.R
# lists what is left.
list(
  "NEWS.md" = list(
    c("  brms's Stan functions where those are accurate. brms's are not
  everywhere: its Frechet log survivor is off by up to 3.8e-3 and is
  `-Inf` in the far right tail, and its GEV density loses digits near
  `xi = 0`, 6.8e-8 at `xi = 1e-7` (`dev/famuni-findings.md`).",
      "  brms's Stan functions where those are accurate. brms's are not
  everywhere. Its Frechet log survivor is off by up to 3.8e-3 and is
  `-Inf` in the far right tail. Its GEV density loses digits near
  `xi = 0`, 6.8e-8 at `xi = 1e-7` (`dev/famuni-findings.md`)."),
    c("    0.001. `cens()` follows frmtmb's inclusive discrete convention, so
    right and interval censoring differ from brms, which reads a right
    bound as `P(Y > y)`; subtract one from those bounds to reproduce a
    brms fit.",
      "    0.001. `cens()` follows frmtmb's inclusive discrete convention.
    So right and interval censoring differ from brms, which reads a
    right bound as `P(Y > y)`. Subtract one from those bounds to
    reproduce a brms fit."),
    c("* `family = \"wiener\"` (or a brms `wiener()` object) takes the family
  from frmtmb.eam when that package is loaded, and otherwise is refused
  with a message that names frmtmb.eam.",
      "* `family = \"wiener\"` (or a brms `wiener()` object) takes the family
  from frmtmb.eam when that package is loaded. Otherwise it is refused
  with a message that names frmtmb.eam.")),
  "extensions/frmtmb.sample/NEWS.md" = list(
    c("* `frm_sample()` samples a `gen_extreme_value()` `xi` without a formula
  as brms samples it, through `tmp_xi`, mapped into the support each
  draw's `mu` and `sigma` leave, with brms's default `normal(0, 2.5)` on
  `tmp_xi`.",
      "* `frm_sample()` samples a `gen_extreme_value()` `xi` without a formula
  as brms samples it, through `tmp_xi`. It maps `tmp_xi` into the
  support that each draw's `mu` and `sigma` leave, and gives it brms's
  default `normal(0, 2.5)`."),
    c("  sampled space. On the review's n = 40 data sets at `xi` -0.45, 0.20
  and -0.85 (4 chains x 2000), 1187, 263 and 3856 of 4000 transitions
  diverged, with Rhat 2.30 on the last, against 0, 1 and 0 in brms.",
      "  sampled space. The review fitted n = 40 data sets at `xi` -0.45,
  0.20 and -0.85 (4 chains x 2000). There 1187, 263 and 3856 of 4000
  transitions diverged, with Rhat 2.30 on the last, against 0, 1 and 0
  in brms."),
    c("  The hold is material where the data put mass there: on the xi -0.85
  data brms puts 29 percent of its posterior below -1, and
  `frm_sample()`'s posterior is brms's conditioned on `xi > -1` (median
  of `xi` -0.876 against -0.877; the review's round 2).",
      "  The hold is material where the data put mass there. On the xi -0.85
  data brms puts 29 percent of its posterior below -1. There
  `frm_sample()`'s posterior is brms's conditioned on `xi > -1` (median
  of `xi` -0.876 against -0.877; the review's round 2).")),
  "extensions/frmtmb.eam/NEWS.md" = list(
    c("  Under an identity link the density, and the fitted objective, is `NaN`
  for a `bs` at or below zero, a `bias` outside (0, 1), or an `ndt` at
  or below zero (`link_ndt = \"identity\"`), where the series alone
  returned a finite number for a bias of 1.2.",
      "  Under an identity link the density and the fitted objective are
  `NaN` outside the support. That is a `bs` at or below zero, a `bias`
  outside (0, 1), or an `ndt` at or below zero
  (`link_ndt = \"identity\"`). Before, the series alone returned a finite
  number for a bias of 1.2."),
    c("  then on brms's scale. The default stays the bounded link. On a scalar
  `ndt` the two links are the same model, and over 200 simulated data
  sets the two maxima agreed to 6e-10 in the log-likelihood.",
      "  then on brms's scale. The default stays the bounded link. On a scalar
  `ndt` the two links are the same model. Over 200 simulated data sets
  the two maxima agreed to 6e-10 in the log-likelihood."))
)
