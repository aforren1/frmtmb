# Query the feature compatibility registry

Gives the declared status of a pair of package features. The status is
one of:

## Usage

``` r
frm_compat(feature_a = NULL, feature_b = NULL, status = NULL)
```

## Arguments

- feature_a, feature_b:

  Feature names, as given by
  [`frm_compat_features()`](https://aforren1.github.io/frmtmb/reference/frm_compat_features.md).
  Supply both for one pair, one for every pair involving that feature,
  or neither for the whole table. Both accept a vector, which gives
  every pair in the cross of the two sides; an empty vector is an error
  rather than an empty answer.

- status:

  Optional character vector; keep only these statuses.

## Value

A data frame with columns `feature_a`, `kind_a`, `feature_b`, `kind_b`,
`status`, and `note`. Family pairs are omitted, because a model carries
one family.

## Details

- `works`:

  The combination is supported and exercised.

- `conditional`:

  Supported, but the note states a condition that the combination must
  meet.

- `refused`:

  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) or the
  post-fit method stops with an error. The refusal is deliberate.

- `broken`:

  The combination is accepted but the result is wrong, or it fails with
  an error that does not explain itself. Avoid the pair. The note gives
  the evidence.

- `untested`:

  Nothing checks this pair. It may work. Treat a silent success as
  unverified, not as support.

The last status is the point of the registry. A guard that does not
exist looks exactly like a guard that passed, so the absence of an error
was never evidence of support.

## See also

[`frm_compat_rules()`](https://aforren1.github.io/frmtmb/reference/frm_compat_rules.md),
[`frm_compat_features()`](https://aforren1.github.io/frmtmb/reference/frm_compat_features.md)

## Examples

``` r
# one pair
frm_compat("rescor", "cens()")
#>   feature_a    kind_a feature_b kind_b  status
#> 1    rescor structure    cens()  aterm refused
#>                                                                        note
#> 1 Refused. This pair was once accepted with the censoring silently dropped.

# everything known about truncation
frm_compat("trunc()", status = c("refused", "broken"))
#>    feature_a kind_a                 feature_b    kind_b  status
#> 1    trunc()  aterm                   student    family refused
#> 2    trunc()  aterm         shifted_lognormal    family refused
#> 3    trunc()  aterm               skew_normal    family refused
#> 4    trunc()  aterm                exgaussian    family refused
#> 5    trunc()  aterm              asym_laplace    family refused
#> 6    trunc()  aterm                     Gamma    family refused
#> 7    trunc()  aterm                      beta    family refused
#> 8    trunc()  aterm                   tweedie    family refused
#> 9    trunc()  aterm               negbinomial    family refused
#> 10   trunc()  aterm                   nbinom1    family refused
#> 11   trunc()  aterm                 geometric    family refused
#> 12   trunc()  aterm                   compois    family refused
#> 13   trunc()  aterm                  binomial    family refused
#> 14   trunc()  aterm                 bernoulli    family refused
#> 15   trunc()  aterm             beta_binomial    family refused
#> 16   trunc()  aterm               multinomial    family refused
#> 17   trunc()  aterm     zero_inflated_poisson    family refused
#> 18   trunc()  aterm zero_inflated_negbinomial    family refused
#> 19   trunc()  aterm    zero_inflated_binomial    family refused
#> 20   trunc()  aterm        zero_inflated_beta    family refused
#> 21   trunc()  aterm            hurdle_poisson    family refused
#> 22   trunc()  aterm              hurdle_gamma    family refused
#> 23   trunc()  aterm          hurdle_lognormal    family refused
#> 24   trunc()  aterm                cumulative    family refused
#> 25   trunc()  aterm                    sratio    family refused
#> 26   trunc()  aterm                    cratio    family refused
#> 27   trunc()  aterm                      acat    family refused
#> 28   trunc()  aterm               categorical    family refused
#> 29   trunc()  aterm                 von_mises    family refused
#> 30   trunc()  aterm                      mi()     aterm refused
#> 31   trunc()  aterm                      ar()   autocor refused
#> 32   trunc()  aterm                      ma()   autocor refused
#> 33   trunc()  aterm                    arma()   autocor refused
#> 34   trunc()  aterm                    cosy()   autocor refused
#> 35   trunc()  aterm                   unstr()   autocor refused
#> 36   trunc()  aterm                quadrature      mode refused
#> 37   trunc()  aterm                    rescor structure refused
#> 38   trunc()  aterm                   mixture structure refused
#> 39   trunc()  aterm               mixture_mvn structure refused
#> 40   trunc()  aterm              bar_crossing   grammar refused
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               note
#> 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Refused: trunc() needs a family with an AD log-CDF.
#> 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Refused: trunc() needs a family with an AD log-CDF.
#> 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Refused: trunc() needs a family with an AD log-CDF.
#> 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Refused: trunc() needs a family with an AD log-CDF.
#> 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Refused: trunc() needs a family with an AD log-CDF.
#> 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Refused: trunc() needs a family with an AD log-CDF.
#> 7                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Refused: trunc() needs a family with an AD log-CDF.
#> 8                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Refused: trunc() needs a family with an AD log-CDF.
#> 9                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Refused: trunc() needs a family with an AD log-CDF.
#> 10                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 11                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 12                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 13                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 14                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 15                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 16                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 17                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 18                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 19                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 20                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 21                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 22                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 23                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 24                                                                                                                                                                                                                                                                                                                                                                                                                                                          Refused: ordinal families carry no AD log-CDF over the response scale.
#> 25                                                                                                                                                                                                                                                                                                                                                                                                                                                          Refused: ordinal families carry no AD log-CDF over the response scale.
#> 26                                                                                                                                                                                                                                                                                                                                                                                                                                                          Refused: ordinal families carry no AD log-CDF over the response scale.
#> 27                                                                                                                                                                                                                                                                                                                                                                                                                                                          Refused: ordinal families carry no AD log-CDF over the response scale.
#> 28                                                                                                                                                                                                                                                                                                                                                                                                                                                            Refused for the same reason: no order, no CDF, no truncation window.
#> 29                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Refused: trunc() needs a family with an AD log-CDF.
#> 30                                                                                                                                                                                                                                                                                                                                                                                                                                            Refused: mi() cannot be combined with cens(), trunc(), or se() on the same response.
#> 31                                                                                                                                                                                Refused: the group's density is joint, so there is no per-row contribution for a frequency weight to repeat, a censoring indicator to replace with a tail probability, a truncation bound to renormalize, or a known standard error to add to. brms refuses weights(), cens() and trunc() here with 'Invalid addition arguments for this model'.
#> 32                                                                                                                                                                                Refused: the group's density is joint, so there is no per-row contribution for a frequency weight to repeat, a censoring indicator to replace with a tail probability, a truncation bound to renormalize, or a known standard error to add to. brms refuses weights(), cens() and trunc() here with 'Invalid addition arguments for this model'.
#> 33                                                                                                                                                                                Refused: the group's density is joint, so there is no per-row contribution for a frequency weight to repeat, a censoring indicator to replace with a tail probability, a truncation bound to renormalize, or a known standard error to add to. brms refuses weights(), cens() and trunc() here with 'Invalid addition arguments for this model'.
#> 34                                                                                                                                                                                Refused: the group's density is joint, so there is no per-row contribution for a frequency weight to repeat, a censoring indicator to replace with a tail probability, a truncation bound to renormalize, or a known standard error to add to. brms refuses weights(), cens() and trunc() here with 'Invalid addition arguments for this model'.
#> 35                                                                                                                                                                                Refused: the group's density is joint, so there is no per-row contribution for a frequency weight to repeat, a censoring indicator to replace with a tail probability, a truncation bound to renormalize, or a known standard error to add to. brms refuses weights(), cens() and trunc() here with 'Invalid addition arguments for this model'.
#> 36 Refused: the truncation normalizer is log(F(ub) - F(lb)) over plain CDFs, and the Gauss-Kronrod nodes reach random-effect values where that difference underflows to exactly zero while the density is still representable. The integrand is then +Inf and the marginalized objective is -Inf, at the Laplace optimum as well as at the starting values, so the fit used to report logLik = +Inf as converged. Laplace stays near the mode and is unaffected: use quadrature = FALSE, REML, or profile for truncated responses.
#> 37                                                                                                                                                                                                                                                                                                                                                                                                                                                      Refused. This pair was once accepted with the truncation silently dropped.
#> 38                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  Refused: mixture() has no CDF.
#> 39                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Refused: mixture_mvn() has no CDF.
#> 40                                                                                                                                                                                                                                                                                                    Refused: a bar term crossed with * or : (as in x * (1 | g)) is not a random-effect specification (lme4#196). Write the crossing inside the bar: (x | g). This spelling was once accepted with the crossing silently dropped.

# the pairs to avoid
frm_compat(status = "broken")[, 1:5]
#> [1] feature_a kind_a    feature_b kind_b    status   
#> <0 rows> (or 0-length row.names)
```
