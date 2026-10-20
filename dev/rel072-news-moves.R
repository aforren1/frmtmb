# The release review (dev/reviews/2026-10-20-release.md, finding 3 and
# section 7): four behavior changes the lanes filed under "New
# features" or "Bug fixes" belong under "Breaking changes". Each entry
# removes a fixed string from the top section (matched once) and adds
# a bullet at the end of that section's "Breaking changes".
# dev/rel072-news-merge.R sources this after the splits
# (dev/rel072-news-split.R), so the strings are the split text.
list(
  "NEWS.md" = list(
    list(remove = r"--(A brms `categorical(refcat = "b")` object passed as
  `family =` keeps its `refcat`; 0.71.0 dropped it and fitted reference
  `a` under the names `mub`, `muc`. )--",
         bullet = r"--(* **A brms `categorical(refcat = "b")` object passed as `family =`
  keeps its `refcat`.** 0.71.0 dropped it and fitted reference `a`
  under the names `mub` and `muc`. The same call now fits brms's
  parameterization, with brms's names.)--"),
    list(remove = r"--(* `validate_prior()` and `prior_summary()` list the class `"sd"` rows in
  brms's order: the class-wide rows, then group by group, with the
  responses inside each group.

)--",
         bullet = r"--(* **`validate_prior()` and `prior_summary()` list the class `"sd"`
  rows in brms's order**: the class-wide rows, then group by group,
  with the responses inside each group.)--"),
    list(remove = r"--(
  `conditional_effects(categorical = TRUE)` on `categorical()` draws its
  Wald band from `fitted()`, where it asked for `band = "boot"`.)--",
         bullet = r"--(* **`conditional_effects(categorical = TRUE)` on `categorical()` draws
  a Wald band from `fitted()` by default**, where it asked for
  `band = "boot"`.)--"),
    list(remove = r"--(* When the standard-error check removes nothing, a log sd with a small
  curvature is interior. The covariance is then taken in the design's
  coordinates. Before, `solve()` of the natural Hessian kept it
  1.4e-3 from the exact one at degree 7.

)--",
         bullet = r"--(* **When the standard-error check removes nothing, a log sd with a
  small curvature is interior.** The covariance is then taken in the
  design's coordinates, which moves the standard errors of such fits.
  Before, `solve()` of the natural Hessian kept it 1.4e-3 from the
  exact one at degree 7.)--"))
)
