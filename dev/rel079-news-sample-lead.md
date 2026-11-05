Needs frmtmb 0.79.0. The exact `gp()` field is non-centered in core, so the
sampler has nothing to transform. It reads the sampling-API exports
`gp_tape_scope()`, `gp_tape_order()`, `gp_nugget_message()`, and the new
readers of a `me()` draw and of a latent residual.
