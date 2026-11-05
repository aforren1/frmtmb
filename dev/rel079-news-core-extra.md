
## New features

* `?"frmtmb-gp"` is a new help topic. It states the nugget, the warning of
  class `frmtmb_gp_nugget`, what to do about it, and the order of the field.

* `frmtmb_control(check_mi = c("warning", "ignore", "stop"))` sets what
  `frm()` does when the Laplace approximation of a read `mi()` value moves
  an estimate.

* `options(frmtmb.latres_dense_rows = )` sets the row limit above which a
  standard error refuses to read the whole joint covariance of a latent
  residual. The limit is 20000 rows. The size is `16 n^2` bytes: about 200 GB
  at 100000 rows, where `conditional_effects()` used 3 GB at 12000 rows.

## Bug fixes

* `gp(x, k = )` has no nugget, but it is not a remedy by itself: `k` must be
  large enough for the field. On `sin(8 * x)` over 0 to 6 at noise 1e-3, `k =
  30` stops at frequency 6.5, reads `sigma` 720 times high and reports
  convergence, with no warning; `k = 60` reads brms's `sigma`.

* A gap of 94 to 1127 nats between `frm()` and brms on an exact `gp()`, read
  earlier, was in the test harness: frmtmb sorts a term's positions and brms
  keeps first appearance, and the translator did not reorder them.

* `fitted()` on a new row of a latent residual holds the latent's draw in
  `Est.Error`, and its `Q` columns are read from that draw, as on a `me()` row.
  On poisson `Est.Error` read 0.133 where the draws' standard deviation is
  1.005. A row with both a `me()` value and a latent residual reads one mean
  and one variance over the two.

* `(1 | p | gr(g, cor = FALSE))` beside `(1 | p | g)` builds two uncorrelated
  blocks, as brms does, in either order. `frm()` refused it. A key written `id
  =` is still refused, since brms reads the first term's flag there.

* `frm_curve()` at its default reads a row of `NA` as `NA`, and the other rows
  answer. Before, it stopped with "infinite or missing values in 'x'".
