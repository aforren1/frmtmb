# Lane parity2: add the CLOSED lines under the items the lane closed and
# the lane's filed items, to dev/test-backlog.md (CRLF kept). Each CLOSED
# line goes directly under its item's last line (line numbers at
# ec8178d3), inserted bottom-up so earlier numbers do not move. Run once.
$p = "C:/Users/adf44/source/r/frmtmb-wt-parity2/dev/test-backlog.md"
$lines = New-Object System.Collections.Generic.List[string]
$lines.AddRange([string[]](Get-Content -LiteralPath $p -Encoding UTF8))
$pre = "  CLOSED in the round of 2026-11-13 (lane parity2): "
$closed = [ordered]@{
  5346 = "the term keeps the unused level (N_1 = 7, as brms); the zero fixed-effect column stays dropped, filed as a decision."
  5336 = "loo_epred(), loo_linpred(), loo_predict(), loo_predictive_interval() and loo_R2() take a subset of the fitted rows."
  5328 = "an ordinal posterior_epred() names its draw margin ""1""..""S"" as brms; categorical stays NULL."
  5324 = "get_predicted() refuses a multinomial or simplex fit by name (brms fails inside insight)."
  5318 = "a 0 + Intercept or center = FALSE intercept keeps its place among the b_, as brms orders it."
  5305 = "frm_sample() draws carry err[k] (sderr * zerr under cov = FALSE, the residual otherwise); the corrected estimate, get_joint_cov() on a new row and M7 are filed again."
  4340 = "default_prior() lists an rr() model (no sd row for the block)."
  4162 = "attr(x, ""min"") of an integer mo() variable sets the lowest category at the fit (brms's Jmo = 4 on 1..4 with min 0); new data are coded from it."
  3924 = "default effects follow brms's get_all_effects() over every dpar and nlpar; names equal brms's on 11 formulas (cs() only in the per-category key)."
  3921 = "s(exp(x)) fits (bs_sexpx_1, sds_sexpx_1); the frame evaluates a smooth's term as an expression."
  3375 = "the REML df counts the fixed effects, as lme4's; AIC 1755.628 on sleepstudy as lme4's."
  2981 = "the j > 1 imputations take the first one's grid in place of the caller's int_conditions."
  1563 = "a later lf() or nlf() replaces the parameter with brms's message; update(fit, bf()) pooling is filed again."
  1502 = "band = ""boot"" and draws lend an observed level's effects by default (sample_new_levels = ""uncertainty""); ""gaussian"" and ""old_levels"" are taken."
}
foreach ($e in ($closed.GetEnumerator() | Sort-Object { [int]$_.Key } -Descending)) {
  $lines.Insert([int]$e.Key, $pre + $e.Value)
}
$i = $lines.IndexOf("## Filed at the 0.82.0 release (2026-11-11)")
if ($i -lt 0) { throw "anchor not found" }
$sec = $lines.IndexOf("## Filed in the round of 2026-11-13")
$new = @(
  "### From lane parity2",
  "",
  "- **LOW. ``update(fit, bf(y ~ x2))`` with a complete ``bf()`` replaces the",
  "  model**, where brms pools the old parameter formulas into it (the half",
  "  of 1560 not done in the round of 2026-11-13).",
  "- **LOW. ``bf(y ~ x, sigma ~ z, sigma ~ w)`` and ``bf(f, sigma ~ w)`` on an",
  "  ``f`` that sets ``sigma`` are refused** (""Duplicated dpar formula""), where",
  "  brms 2.23.0's ``brmsformula()`` keeps the later one with ""Replacing",
  "  initial definitions of parameters"". A later ``+ lf()`` does that since",
  "  this round (``bform_replace_pars()``); ``bf_dots()`` does not.",
  "- **LOW. The rest of the latent residual bundle (5087, 5300)**: a",
  "  corrected estimate and logLik (Laplace + C, ``theta + shift``),",
  "  ``get_joint_cov()`` on a new row's design, which reads ``V`` of the outer",
  "  parameters alone, and M7. Only ``err[k]`` was done (dev/parity2-findings.md).",
  "- **LOW (decision). ``frm(drop_unused_levels = FALSE)`` drops the zero",
  "  column of an unused factor level** (``fc`` of ``y ~ f``), where brms keeps it",
  "  (K = 3, ``b_fc`` drawn from its prior, improper under brms's flat default).",
  "  A coefficient of a zero column has no maximum likelihood estimate.",
  "  Recommendation: keep the drop. dev/parity2-findings.md, item 5337.",
  "- **LOW. ``conditional_effects()`` leaves a factor-smooth's own grouping",
  "  factor out of its default displays** (``s(x, g, bs = ""fs"")``), where",
  "  brms's ``get_all_effects()`` lists ``g`` and ``x:g``. frmtmb drops that term",
  "  at the population level, so the display would draw one line per level",
  "  of a term the curve does not read. Kept; recorded for the ledger.",
  "- **LOW. The default displays of ``mo()``, ``me()`` and ``mi()`` terms come",
  "  in that order**, where brms lists its special terms in written order. The",
  "  set is brms's; only the order of the list can differ when a formula",
  "  mixes them.",
  ""
)
if ($sec -lt 0) {
  $hdr = @("## Filed in the round of 2026-11-13", "")
  $lines.InsertRange($i, [string[]]($hdr + $new))
} else {
  $lines.InsertRange($i, [string[]]$new)
}
$enc = New-Object System.Text.UTF8Encoding($false)
[System.IO.File]::WriteAllText($p, (($lines -join "`r`n") + "`r`n"), $enc)
"done: " + $lines.Count
