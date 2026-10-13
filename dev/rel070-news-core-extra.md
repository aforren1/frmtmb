## Extension API

* **The sampling API exports `simplex_frame_terms()`**
  (`?frmtmb-sampling-api`): every simplex of a fit's frame in the
  shape `mo_frame_terms()` gives a `mo()` term, the `mo()` simplexes
  and then each `cox()` baseline under brms's name `sbhaz[k]`.
  `mo_chart_frame()` now sets a `cox()` baseline's chart too.
  frmtmb.sample 0.18.0 samples every simplex through them.
