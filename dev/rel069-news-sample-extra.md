## Bug fixes

* `check_laplace()` leaves a parameter whose `z_shift` or `sd_ratio` is
  not finite out of its "Laplace/Wald approximation questionable"
  message. The test passed `NA` to `if (any(...))`, so such a parameter
  stopped the call ("missing value where TRUE/FALSE needed") when no
  other parameter was flagged, and was listed as `NA` when one was. The
  returned table is unchanged.
