### brms-28. `logistic_normal()` in `mvbf()` does not compile

- Seen: 2.23.0, 2026-10-14 (lane catmulti, brms-C;
  `dev/catmulti-brms-fits.R`, case "mv-ln"; the Stan code in its log
  `brms-probe.txt`).
- Repro: `brm(bf(Y ~ x, family = logistic_normal()) + bf(w ~ x, family
  = gaussian()) + set_rescor(FALSE), data = dd)`, `Y` a three-column
  composition.
- Expected: a fit. Actual: "Semantic error ... Identifier "Llncor" not
  in scope. Did you mean "Lncor"?". The generated quantities read
  `Llncor` and `nlncor`, where the parameters and the transformed data
  declare `Llncor_Y` and `nlncor_Y` (`stan-predictor.R`: the glue of
  `gen_def` has no `{p}` suffix). With two `logistic_normal()`
  responses the generated quantities also declare `Lncor` and `lncor`
  twice, and `rename_family_cor_pars()` would give both responses'
  correlations the one name `lncor__b__c`.
- Workaround: frmtmb fits the model. Its test repairs the two
  identifiers and checks brms's model block (lines 135 and 136 of 145
  change).
- Record: `dev/catmulti-findings.md`, sections 8 and 13.
- Status: not reported.

### brms-29. `update()` without new data computes a data-dependent basis again on fewer rows

- Seen: 2.23.0, 2026-10-19 (lane drop1raw, review finding 2;
  `dev/drop1raw-rev-brms2.R`, log `brms2.txt`, run on an empty fit with
  a stub stanfit, so `update()` runs without a compile).
- Repro: a fit of `z ~ poly(x2, 2) + x3` on 40 rows (seed 5), 2 of
  them `NA` in `x3`. The fit's `poly()` column is computed on all 40
  rows before `na.action`; `object$data` keeps the 38 rows after it, and
  `update(fit)` computes `poly()` again on those 38. The script prints
  "fit column == update column: FALSE".
- Expected: an update with no change reproduces the fit's design.
  Actual: a different basis, so the refit is a different model (in
  frmtmb's measurement of the same defect, `ns(x2, 3) + x3` moved the
  log-likelihood by 1.5e-2).
- Workaround: frmtmb keeps the fit's bases (`predvar_map`) on every
  refit from the fit's own data (`update()` without `newdata`,
  `drop1()`, `influence()`), and so does not follow brms here.
- Record: `dev/reviews/2026-10-19-drop1raw.md`, finding 2;
  `dev/drop1raw-findings.md`, fix round 1.
- Status: not reported.
