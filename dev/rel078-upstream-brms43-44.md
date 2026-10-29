### brms-43. One response of a student `set_rescor(TRUE)` model reads `nu = 1`

- Seen: 2.23.0, 2026-10-28 (`dev/rescorpp-brmsstudent.R`, log
  `dev/rescorpp-brmsstudent.log`; found by lane rescorpp,
  `dev/rescorpp-findings.md`, section 2; confirmed by its review,
  `dev/reviews/2026-10-28-rescorpp.md`, "The shared-nu divergence is
  justified").
- Repro: data set `st` of `dev/rescorpp-data.R` (two student responses,
  nu 5, seed 13, 200 rows), fitted by `dev/rescorpp-brmsfit.R` (2
  chains, 2000 iterations, seed 21) as
  `mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE)` with
  `family = student()`.
  Then `log_lik(fit, resp = "y1")` and
  `posterior_predict(fit, resp = "y1")`.
- Expected: one response's density and draw are the margin of the
  multivariate t, a t with the shared `nu`. Actual: the one-response
  preparation has `nu` 1 in every draw (the draws of `nu` are 8.64,
  5.41, ...). `log_lik(resp = "y1")` equals `dt(z, 1)` to a maximum
  absolute difference of 0 and differs from `dt(z, nu)` by up to 6.44.
  `posterior_predict(resp = "y1")` is a Cauchy draw: the absolute
  prediction error at the 50, 90 and 99 percent quantiles is 1.09, 6.95
  and 68.8, against 0.78, 2.15 and 4.19 for y1 of the joint draw.
- Cause: `brmsterms.R` lines 184 to 189 set `pfix$nu <- 1` on every
  response when `no_nu()` is true (rescor with student), a placeholder
  so that the per-response Stan code declares no `nu`.
  `prepare_predictions.mvbrmsframe()` (`prepare_predictions.R` lines 61
  to 95) puts the shared `nu` in `out$dpars$nu` only when
  `length(resp) > 1`; one `resp` takes `x$terms[[resp]]` and so the
  placeholder.
- Workaround: frmtmb.sample 0.25.0 draws and evaluates one response at
  the shared `nu`, a deliberate divergence stated in its NEWS.
- Record: `dev/rescorpp-findings.md`, section 2 and "Filed, not fixed".
  `gh search issues --repo paul-buerkner/brms` found no issue,
  2026-10-28.
- Status: not reported.

### brms-44. A `resp` out of model order stops every method that prepares several responses under `set_rescor(TRUE)`

- Seen: 2.23.0, 2026-10-28 (`dev/rescorpp-compare.R`, log
  `dev/rescorpp-compare-g3.log`; the review's
  `dev/review-rescorpp-layout.R`, logs
  `dev/review-rescorpp-layout-*.log`; lane rescorpp,
  `dev/rescorpp-findings.md`, section 2, corrected after the review's
  finding 8).
- Repro: data set `g3` of `dev/rescorpp-data.R` (three gaussian
  responses, seed 11, 200 rows), fitted by `dev/rescorpp-brmsfit.R` as
  `mvbf(bf(y1 ~ x), bf(y2 ~ x), bf(y3 ~ x)) + set_rescor(TRUE)`. Then
  any of
  `posterior_predict()`, `posterior_epred()`, `posterior_linpred()`,
  `predictive_error()`, `fitted()`, `predict()`, `residuals()`,
  `log_lik()` and `loo_predict()` with `resp = c("y2", "y1")` (or
  `c("y3", "y1")`).
- Expected: the responses named, in the order named, with the
  correlations of that order, or the model's order with a message.
  Actual: every call stops with "The following variables are missing in
  the draws object: {'rescor__y2__y1'}", `posterior_epred()` too, which
  reads no residual correlation. The same on the `gs`, `st` and `re`
  fits of the review. Under `set_rescor(FALSE)` the calls answer.
- Cause: `get_cornames(resp, ...)` in
  `prepare_predictions.mvbrmsframe()` builds the correlation names in
  the order given, and the draws hold them in the model's order only.
- Workaround: frmtmb.sample 0.25.0 answers a permuted `resp`, with the
  correlations in the order named, a divergence stated in its NEWS.
- Record: `dev/rescorpp-findings.md`, section 2 and "Filed, not fixed";
  `dev/reviews/2026-10-28-rescorpp.md`, finding 8. No upstream issue
  found by `gh search issues` on `rescor`, 2026-10-28.
- Status: not reported.
