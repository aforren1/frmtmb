# Contribute to the compatibility matrix from another package

[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
answers what one feature does in the presence of another, and it answers
in three states: the pair works, the pair is refused, or the pair is
untested. A feature that lives outside frmtmb has to be able to say the
same things about itself, or the matrix reports on part of the package
and reads as if it reported on all of it.

## Usage

``` r
frmtmb_register_compat(features = NULL, rules = NULL, expects = character(0))

compat_rule_builder()
```

## Arguments

- features:

  Named character vector mapping a feature's DISPLAY name to its kind:
  `c("hmm" = "structure")`, `c("dec()" = "aterm")`. The display name is
  how the feature is written in a formula or a call, parentheses
  included for a callable one; the kind groups it in the printed matrix
  and decides which kind-level rules reach it, and is one of `"family"`,
  `"covstruct"`, `"aterm"`, `"special"`, `"autocor"`, `"mode"`,
  `"structure"`, `"method"`, `"grammar"`. A name the registry already
  carries under the same kind is a no-op, so a reload is not an error
  and neither is declaring an addition term that
  [`frmtmb_register_aterm()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_aterm.md)
  has already declared for you.

- rules:

  A FUNCTION of no arguments returning a rule data frame, not the data
  frame itself, so that contributed rules are built on demand exactly as
  the core ones are. Build the frame with `compat_rule_builder()`: it
  returns a list of two functions,
  `r(a, b, status, note, override = FALSE)` to record one rule and
  `rules()` to return the accumulated frame. `a` and `b` are display
  names, `"*"` matches every feature, and `override = TRUE` says the
  rule beats a more specific one rather than losing to it.

- expects:

  Character vector of feature names the rules refer to but this package
  does not supply, because another package does. They are exempt from
  the check below and are NOT added to the vocabulary: a rule naming one
  of them lies dormant until the package that owns the feature is
  loaded, which is what makes the rule true. Use it for that and nothing
  else. A feature in the table that nothing implements is the failure
  the registry exists to prevent.

  A name the vocabulary already holds is refused, because the
  declaration then exempts nothing: give it a kind in `features =`, or
  drop it. A name still unresolved when
  [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
  is called is reported there, in the `unresolved` attribute the print
  method shows, so a rule waiting for a package that never arrives is
  visible rather than silent. That is what makes a misspelling inside
  `expects` findable: the check here cannot catch one, since the whole
  point is that the owner is absent.

  The near-miss guess in a refusal is searched over the vocabulary only,
  never over `expects`, so a misspelling of an expected name is refused
  without a suggestion. Pointing at a feature that is not in the session
  would be its own confusion.

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

## What registration refuses

A rule side names a feature, a whole kind (`"kind:family"`), a named
group (`"group:cdf"`), or everything (`"*"`). A side that names none of
those matches no pair, and a rule that matches no pair is dropped
without a word: the contributing package's table then reads as complete
while covering less than it claims, which is the worst thing a
compatibility surface can do. So registration refuses it, naming the
rule, the spelling, and the nearest entry in the vocabulary when one is
close.

The check runs against the vocabulary AFTER `features` is added, so a
package names its own new features freely in its own rules. It runs at
registration, which means `rules` is called once there: keep it a pure
builder that reads nothing but its own arguments.

An addition term registered with
[`frmtmb_register_aterm()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_aterm.md)
is added to the vocabulary for you, as `"<name>()"` of kind `"aterm"`. A
registered term the table cannot describe would be a gap by
construction, so the registrant is not asked to say it twice; saying it
anyway is a no-op. Register the term BEFORE the rules that name it.

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
  b$r("wiener", "hmm", "untested",
      "Another package's feature, so the rule waits for it.")
  frmtmb_register_compat(features = c("wiener" = "family"),
                         rules = b$rules, expects = "hmm")
}
# the accumulator on its own, which is all a rule set is
b <- compat_rule_builder()
b$r("dec()", "trials()", "refused", "Different response shapes.")
b$rules()
#>   feature_a feature_b  status                       note override
#> 1     dec()  trials() refused Different response shapes.    FALSE
```
