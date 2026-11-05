# prepends the development sections; run once
pre <- function(file, txt) {
  x <- readChar(file, file.size(file), useBytes = TRUE)
  if (grepl("\r\n", x, fixed = TRUE, useBytes = TRUE)) txt <- gsub("\n", "\r\n", txt, fixed = TRUE)
  writeBin(charToRaw(paste0(txt, x)), file)
}
pre("NEWS.md", "# frmtmb (development version)

## Breaking changes

* The exact `gp()` is non-centered: `b` of the block is z, a standard normal,
  brms's `zgp`, and the field is `sdgp * L(lscale) * z`, with `L` the Cholesky
  factor of the kernel plus a nugget of `1e-12` on the correlation. The field
  was centered at a nugget of `1e-6`, whose white noise of variance
  `1e-6 * sdgp^2` came off `sigma^2`, or moved `sdgp` and the length scale when
  positions repeat, `sigma` is known or the family has none. Every fit that
  relied on it moves. On a curve `x^2` over x from 0 to 10 with `sigma` 0.3,
  `frm()` read `sigma` 6e-5 and now reads 0.35, the marginal's maximum. On
  `sin(2 * x)` with noise 1e-3 it read `sigma` 1.3e-7 and now reads 1.1e-3, 1.8
  nats higher. At a known `se(3e-4)` it was 29 nats short. On `exp(0.8 * x)`
  with `sigma` 0.3 the fit was 72 nats short and is now 0.7 short, with a
  warning that it stopped before a maximum: the kernel's condition there
  (2e12 against the noise) leaves the objective rounding at 1e-3.
* The nugget is relative to `sdgp`, so the fit of `y * c` is the fit of `y`
  scaled, on any `c`. brms's nugget is an absolute `1e-12`: at `y * 1000` its
  kernel does not factor at the estimates and its chains do not mix (Rhat 3.77),
  and at `y / 1000` it floors `sigma`. At `sdgp` near 1 the two are one number.
* `fit$estimates$b`, `fit$obj`, the joint precision and `start = list(b = )`
  read z in an exact `gp()` block, in sorted position order. `ranef()`, `coef()`
  and `getME(fit, \"b\")` still read the field. Read the field off those.
* An exact `gp()` starts its length scale at the best of four values read from
  the objective, from the largest distance between two positions down to 1/64 of
  it, unless `start = list(theta = )` is given. So the start does not depend on
  x's units, and a field that is rough against its positions finds its maximum:
  `sin(8 * x)` over 0 to 10 with `sigma` 0.3 stopped 43.7 nats short on 1 seed
  of 2, at the maximum with an empty field.
* `frm_curve_deriv()` on a fit with an exact `gp()` takes a default step of
  `5e-4` of the grid's range at order 2, where it was `1e-4`. At the smaller
  nugget a second difference meets the rounding of the kriging row: at `1e-4` a
  standard error read up to 20 percent off.

## Costs

* A fit with an exact `gp()` costs more than the centered one did: measured
  times are in the help of `gp()`'s term under `?frm`. The inner Hessian is
  `I + sd^2 L' W L`, dense in the positions, and TMB's tape reads it column by
  column, so `gr` grows faster than `fn`.

## New features

* `gp_tape_scope()`, in the sampling API: the tape settings a frame with an exact
  `gp()` needs while its tape is built. They change no value.

")
pre("extensions/frmtmb.sample/NEWS.md", "# frmtmb.sample (development version)

## Breaking changes

* An exact `gp()` is non-centered in core, as brms's `zgp`, at a nugget of
  `1e-12 * sdgp^2`, so `frm_sample()` samples the tape `frm()` fits, with the
  sampler's priors added, and has nothing to transform. The block no longer has
  a factor for `reparameterize =`, and that argument does not reach it: an
  exact field is non-centered with or without it. The `b[i]` columns of its
  draws are z in sorted position order, where they were the field. `ranef()`
  and every predictive method read the field through `expand_b()`, as before.

")
