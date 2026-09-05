# Compatibility rules, before resolution

The registry stores rules, not one row per pair. A rule side is one of
four patterns, from least to most specific: `"*"` for any feature
(specificity 0), `"kind:<kind>"` for a whole kind (1), `"group:<group>"`
for a named feature set (2), or a bare feature name (3).

## Usage

``` r
frm_compat_rules()
```

## Value

A data frame with columns `feature_a`, `feature_b`, `status`, `note`,
and `override`.

## Details

A pair takes the status of the rule with the highest precedence, where
precedence is the two sides' specificities sorted descending and
compared lexicographically: `(3,3)` beats `(3,2)` beats `(3,1)` beats
`(3,0)` beats `(2,2)` beats `(2,1)`, and so on. A rule that is strictly
more specific on one side and no less specific on the other therefore
always wins. Rules can tie only when their two specificities match
exactly; the later rule wins those, so an override is appended rather
than inserted, and a tie between rules that disagree about the status is
a registry defect the test suite rejects.

Use
[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
to read the resolved answer for a pair. Use this function to see which
rule is doing the work.

## See also

[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md),
[`frm_compat_features()`](https://aforren1.github.io/frmtmb/reference/frm_compat_features.md)

## Examples

``` r
rules <- frm_compat_rules()
nrow(rules)
#> [1] 325
subset(rules, status == "broken")[, c("feature_a", "feature_b")]
#> [1] feature_a feature_b
#> <0 rows> (or 0-length row.names)
```
