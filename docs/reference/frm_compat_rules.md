# Compatibility rules, before resolution

The registry stores rules, not one row per pair. A rule side is one of
four patterns, from least to most specific: `"*"` for any feature,
`"kind:<kind>"` for a whole kind, `"group:<group>"` for a named feature
set, or a bare feature name. A pair takes the status of the most
specific matching rule; ties go to the later rule.

## Usage

``` r
frm_compat_rules()
```

## Value

A data frame with columns `feature_a`, `feature_b`, `status`, and
`note`.

## Details

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
#> [1] 228
subset(rules, status == "broken")[, c("feature_a", "feature_b")]
#> [1] feature_a feature_b
#> <0 rows> (or 0-length row.names)
```
