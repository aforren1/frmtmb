# Contribute to the compatibility matrix from another package

[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
answers what one feature does in the presence of another, and it answers
in three states: the pair works, the pair is refused, or the pair is
untested. A feature that lives outside frmtmb has to be able to say the
same things about itself, or the matrix reports on part of the package
and reads as if it reported on all of it.

## Usage

``` r
frmtmb_register_compat(features = NULL, rules = NULL)

compat_rule_builder()
```

## Arguments

- features:

  Named character vector mapping a feature's DISPLAY name to its kind:
  `c("hmm" = "structure")`, `c("dec()" = "aterm")`. The display name is
  how the feature is written in a formula or a call, parentheses
  included for a callable one; the kind groups it in the printed matrix
  (`"family"`, `"covstruct"`, `"aterm"`, `"special"`, `"autocor"`,
  `"mode"`, `"structure"`).

- rules:

  A FUNCTION of no arguments returning a rule data frame, not the data
  frame itself, so that contributed rules are built on demand exactly as
  the core ones are. Build the frame with `compat_rule_builder()`: it
  returns a list of two functions,
  `r(a, b, status, note, override = FALSE)` to record one rule and
  `rules()` to return the accumulated frame. `a` and `b` are display
  names, `"*"` matches every feature, and `override = TRUE` says the
  rule beats a more specific one rather than losing to it.

## Value

`NULL`, invisibly. Called for the registration. `compat_rule_builder()`
returns a list with elements `r` and `rules`.

## Details

`frmtmb_register_compat()` is that seam. Call it from the contributing
package's `.onLoad()`: by then every namespace is sealed, so the
collation-order question a top-level call would raise does not arise.
`compat_rule_builder()` is the accumulator the core's own rules are
written with, exported so that a contributed rule reads the same as a
core one.

Contributions are APPENDED, which is the only safe direction. Rules of
equal specificity resolve later-wins, so a contributed rule may override
a core default and a core default can never silently override a
contributed one.

## Status vocabulary

Every rule declares one of three states, and the third is the reason the
registry exists: a guard that does not exist looks exactly like a guard
that passed, so "no error" was never evidence of support.

- `"works"`:

  The pair is supported and exercised by a test.

- `"refused"`:

  The pair is rejected with a message. Say in the note what the user
  should do instead.

- `"untested"`:

  Nothing is known. Neither a promise nor a refusal.

## See also

[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
for the matrix these fill,
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
and
[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
for the family-side seams a contributor usually registers alongside
these, and
[`vignette("compatibility")`](https://aforren1.github.io/frmtmb/articles/compatibility.md)

## Examples

``` r
# what a contributing package's .onLoad() does
contribute <- function() {
  b <- compat_rule_builder()
  b$r("wiener", "cens()", "refused",
      "The family supplies no lcdf, so there is no CDF to censor with.")
  b$r("wiener", "*", "untested",
      "Not exercised outside this package's own suite.")
  frmtmb_register_compat(features = c("wiener" = "family"),
                         rules = b$rules)
}
# the accumulator on its own, which is all a rule set is
b <- compat_rule_builder()
b$r("dec()", "trials()", "refused", "Different response shapes.")
b$rules()
#>   feature_a feature_b  status                       note override
#> 1     dec()  trials() refused Different response shapes.    FALSE
```
