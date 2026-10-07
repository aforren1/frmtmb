## Extension API

* **The sampling API exports `mo_simplex()`, `mo_coords()`,
  `mo_chart_frame()`, `mo_frame_terms()` and `summary_mo_frame()`**
  (`?frmtmb-sampling-api`), so that an extension can read a `mo()`
  simplex in the coordinates the fit holds it in, or in the softmax
  chart a sampler needs, and give the weights brms's `simo_` names.
  frmtmb.sample 0.17.0 samples a `mo()` simplex through them.
