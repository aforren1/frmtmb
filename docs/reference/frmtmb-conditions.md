# Classed conditions

Every error, warning and message that frmtmb raises itself carries a
class that names frmtmb, so a caller can catch all of them without
matching message text:

|           |                                               |
|-----------|-----------------------------------------------|
| raised as | class vector                                  |
| error     | `c("frmtmb_error", "error", "condition")`     |
| warning   | `c("frmtmb_warning", "warning", "condition")` |
| message   | `c("frmtmb_message", "message", "condition")` |

A condition that an extension package raises has one more class in
front, made from the package name with each dot changed to an
underscore. For example, an error from frmtmb.eam has the class
`c("frmtmb_eam_error", "frmtmb_error", "error", "condition")`. Thus
`tryCatch(frmtmb_error = )` catches a refusal from frmtmb and from every
extension, and `tryCatch(frmtmb_eam_error = )` catches the refusals that
are frmtmb.eam's by the rules below.

The subclass names one package, found in this order:

1.  The package that the raising function gives as `package =`. frmtmb
    gives it in two cases. The first is a refusal whose text an
    extension gave as data, such as the `refusals` of
    [`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
    or the `sim_refusal` of
    [`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md).
    The second is a refusal about a declaration of a family: an addition
    term that the family needs or does not read, or a simulator, mean,
    variance, unit deviance, latent state or likelihood factor that the
    family does not supply. In both cases the package is the one that
    built the family or the structure. frmtmb.eam gives it for a refusal
    of its non-decision-time seam
    ([`frmtmb.eam::ndt_bound()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound.html)
    and the functions near it), where it is the package that called the
    seam. Thus a refusal about the bound of
    [`frmtmb.learn::rlddm()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.html)
    is a `frmtmb_learn_error`.

2.  Otherwise, the package of the function that raises the condition.

A family or a structure belongs to the package whose function called
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
or
[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
to build it. No subclass is added for frmtmb itself or for code outside
a package. So a refusal that frmtmb raises for any other reason has no
subclass, also when it is about a model with a family from an extension.
Examples are a refusal of an argument of
[`predict()`](https://rdrr.io/r/stats/predict.html), of a prior, of
`newdata`, or of a value outside the range of a link function. An
extension function that the user calls, such as
[`frmtmb.eam::ndt_time()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_time.html)
on a fit of another package, raises with its own subclass.

This is the contract that brms keeps with its class `brms_error`. frmtmb
builds the conditions with base R, so the class vector does not include
`rlang_error`.

Some conditions do not have these classes:

- An error from another package, such as RTMB, TMB, Matrix or base R,
  keeps its own class. The exception is an error that occurs while
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md)
  optimizes a model: frmtmb catches it and raises it again as a
  `frmtmb_fit_error`, which is also a `frmtmb_error`.

- An argument that the function does not have, such as `frm(priors = )`,
  is refused by R itself with "unused argument".

- The refusal that
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) gives
  for a multivariate fit also has the class `simpleError`, after
  `frmtmb_error`. The insight package examines that class to replace a
  failed standard error with `NA`.

## Usage

``` r
frm_stop(..., call. = TRUE, domain = NULL, class = NULL, package = NULL)

frm_warning(
  ...,
  call. = TRUE,
  immediate. = FALSE,
  noBreaks. = FALSE,
  domain = NULL,
  class = NULL,
  package = NULL
)

frm_message(..., domain = NULL, appendLF = TRUE, class = NULL, package = NULL)

frm_family_package(family)

frm_match_arg(arg, choices, several.ok = FALSE)
```

## Arguments

- ...:

  Zero or more objects that are pasted together, as in
  [`stop()`](https://rdrr.io/r/base/stop.html), or one condition object.

- call.:

  Logical. If `TRUE`, the call of the function that raises the condition
  is part of the condition, as in
  [`stop()`](https://rdrr.io/r/base/stop.html).

- domain:

  Passed to [`.makeMessage()`](https://rdrr.io/r/base/message.html), for
  translation.

- class:

  A character vector of more classes, put in front of the classes that
  frmtmb adds.

- package:

  `NULL`, or the name of the package whose subclass the condition gets.

- immediate.:

  Logical. If `TRUE` and `getOption("warn")` is 0, the warning prints at
  once, as in [`warning()`](https://rdrr.io/r/base/warning.html).

- noBreaks.:

  Ignored. It is accepted so that a call to
  [`warning()`](https://rdrr.io/r/base/warning.html) can use the same
  arguments.

- appendLF:

  Logical. If `TRUE`, a newline is added to the message, as in
  [`message()`](https://rdrr.io/r/base/message.html).

- family:

  A family object from
  [`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
  or
  [`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md).

- arg, choices, several.ok:

  As in [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

## Value

`frm_stop()` does not return. `frm_warning()` returns its message
invisibly. `frm_message()` returns `NULL` invisibly. `frm_match_arg()`
returns the matched choices. `frm_family_package()` returns the name of
the package that built `family`, or `"frmtmb"`.

## For extension authors

`frm_stop()`, `frm_warning()` and `frm_message()` are the functions that
frmtmb uses to raise its conditions. Use them in an extension in place
of [`stop()`](https://rdrr.io/r/base/stop.html),
[`warning()`](https://rdrr.io/r/base/warning.html) and
[`message()`](https://rdrr.io/r/base/message.html). They take the same
arguments as the base functions and give the same message text and the
same call, and they add the classes in the table above. The class of the
extension is found from the namespace of the function that calls them,
so an extension does not supply it. A function of base R between the
two, such as [`lapply()`](https://rdrr.io/r/base/lapply.html) when a
helper is given as its `FUN`, is skipped.

Given one condition object, the helpers signal that object, as the base
functions do, with the frmtmb classes added in front of its own classes.

`class` adds classes in front of the others, for a condition that a
caller must be able to catch by itself. `package` names the package
whose subclass the condition gets, in place of the package of the
calling function. Use it where one package raises a refusal that another
package wrote. `package = "frmtmb"` gives no subclass.
`frm_family_package()` gives the package that built a family, the value
to use for a refusal about that family.

`frm_match_arg()` does what
[`match.arg()`](https://rdrr.io/r/base/match.arg.html) does, and a value
that matches no choice is a `frmtmb_error` that names the argument, the
value and the choices. The subclass is that of the function that calls
`frm_match_arg()`.

These differences from the base functions remain:

- With `call. = TRUE`, an error that nothing catches in a session that
  is not interactive can print one more line, `Calls:`, that includes
  the name of the helper. This line is R's short traceback. The message
  and the call are the same as from
  [`stop()`](https://rdrr.io/r/base/stop.html). frmtmb itself uses
  `call. = FALSE`.

- `noBreaks.` has no effect.

- `frm_message()` records the call to itself, where
  [`message()`](https://rdrr.io/r/base/message.html) records the call to
  [`message()`](https://rdrr.io/r/base/message.html). R does not print
  the call of a message.

## See also

[frmtmb-extension-api](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
for the rest of the interface an extension uses.

## Examples

``` r
d <- data.frame(y = c(1, 2, 3))
e <- tryCatch(frm(y ~ 1, data = d, family = "not_a_family"),
              frmtmb_error = function(e) e)
class(e)
#> [1] "frmtmb_error" "error"        "condition"   
conditionMessage(e)
#> [1] "not_a_family is not a supported family. Supported families are: gaussian, poisson, binomial, Gamma, lognormal, student, negbinomial, nbinom2, nbinom1, beta, Beta, tweedie, compois, zero_inflated_poisson, zero_inflated_negbinomial, hurdle_poisson, multinomial, cumulative, beta_binomial, skew_normal, inverse.gaussian, exgaussian, bernoulli, geometric, exponential, weibull, shifted_lognormal, hurdle_gamma, hurdle_lognormal, zero_inflated_binomial, zero_inflated_beta, asym_laplace, zero_inflated_asym_laplace, huber, sratio, cratio, acat, von_mises, cox, categorical"

f <- function() frm_warning("a warning from f()")
w <- tryCatch(f(), frmtmb_warning = function(w) w)
class(w)
#> [1] "frmtmb_warning" "warning"        "condition"     
conditionCall(w)
#> f()

g <- function(type = c("response", "link")) frm_match_arg(type)
g("resp")
#> [1] "response"
e <- tryCatch(g("bogus"), frmtmb_error = function(e) e)
conditionMessage(e)
#> [1] "`type` must be one of \"response\", \"link\", not character \"bogus\""
```
