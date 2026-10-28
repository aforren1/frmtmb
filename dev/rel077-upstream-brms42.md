### brms-42. `lb` or `ub` on class `cor` gives a Stan program that does not parse

- Seen: 2.23.0, 2026-10-27 (`dev/reslot-brms-corlb.R`, log
  `dev/reslot-out/brms-corlb.txt`; found by lane reslot,
  `dev/reslot-findings.md`, section 1).
- Repro: seed 1, 8 groups of 10, `bf(y ~ x + (1 + x | g))`,
  `stancode(..., prior = set_prior("lkj(2)", class = "cor", lb = 0))`
  or `set_prior("", class = "cor", group = "g", lb = -0.5, ub = 0.5)`,
  then `rstan::stanc(model_code = code)`.
- Expected: `validate_prior()` refuses a bound on class `cor`, since a
  correlation matrix has no scalar bound, or the program parses.
  Actual: `validate_prior()` accepts the bound, and `stancode()` writes
  `cholesky_factor_corr<lower=0>[M_1] L_1` and
  `lkj_corr_cholesky_lccdf(0 | 2)`, which `rstan::stanc()` rejects:
  "Syntax error in 'string', line 46, column 22". The same model
  without the bound parses.
- Workaround: frmtmb takes no bound on a correlation (one prior slot
  per random-effect parameter, 0.77.0), so there is nothing to port.
- Record: `dev/reslot-findings.md`, "Not done, and why".
- Status: not reported.
