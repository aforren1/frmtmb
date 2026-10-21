# The release review's NEWS fixes (dev/reviews/2026-10-21-release.md,
# section 10), as pairs (text after dev/rel073-news-split.R, new text).
# dev/rel073-news-merge.R applies them after the splits, the same way:
# each matched once in the file's top section, any run of whitespace
# standing for a space.
list(
  "NEWS.md" = list(
    # finding 2: on 0.72.0 that model was refused before any fit (no CDF)
    c(r"---[gives no mass to. A `hurdle_negbinomial()` fit with zeros below
      `trunc(lb = 1)` ended at a NaN gradient.]---",
      r"---[gives no mass to.]---"),
    # finding 7: 0.72.0 refused a fractional cens() bound already
    c(r"---[A fractional `cens()` or `trunc()` bound on a count is
      refused.]---",
      r"---[A fractional `trunc()` bound on a count is refused, as a
      fractional `cens()` bound already was.]---"),
    # finding 1, the code fix: the check runs on the uncensored rows,
    # and carries zero_inflated_asym_laplace()
    c(r"---[gen_extreme_value and exgaussian.]---",
      r"---[gen_extreme_value, exgaussian and
      zero_inflated_asym_laplace.]---"),
    c(r"---[The check does not run with `se()`, `cens()`, residual
      correlation, `mi()` or several responses.]---",
      r"---[The check does not run with `se()`, residual correlation,
      `mi()` or several responses. Under `cens()` it counts the
      uncensored rows, because a censored row's probability is at most 1
      and its own effect can bring it to 1. Under the Laplace
      approximation each censored row takes 1/2 off the excess.]---"),
    c(r"---[the effect exactly and is refused only when the weights sum to
      more than the number of rows.]---",
      r"---[the effect exactly and is refused only when the weights sum to
      more than the number of rows. `zero_inflated_asym_laplace()` is
      refused the same way: on its nonzero rows it is asym_laplace
      (log-likelihood 278.05 on 80 unweighted rows). With `cens()` both
      checks count the uncensored rows.]---"),
    # finding 7: the relative clause needs its "that"
    c(r"---[A row rejection has not filled after 20 passes is now drawn by
      inverting its CDF, for every family of this package that has
      one.]---",
      r"---[A row that rejection has not filled after 20 passes is now
      drawn by inverting its CDF, for each family of this package with
      one.]---")),
  "extensions/frmtmb.sample/NEWS.md" = list(
    # finding 3: the regrouping put the six under New features, below
    c(r"---[and those of the six functions above,]---",
      r"---[and those of the six `loo_*()` functions under New
      features,]---"))
)
