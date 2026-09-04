# Make a user-written function tape-safe

Returns `f` with RTMB's tape-safe replacements for
[`c()`](https://rdrr.io/r/base/c.html), `[<-` and `diag<-` bound at the
top of its body. Call it on any function a user supplies that frmtmb
will later put on an AD tape, so that the user does not have to carry
the binding boilerplate. frmtmb applies it to every function slot of
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
and
[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md);
a package that takes tape-side functions of its own applies it to those.
It is safe to call on any object: a non-function and a primitive come
back unchanged.

## Usage

``` r
frmtmb_ad_overload(f)
```

## Arguments

- f:

  A function, or any object. Only a non-primitive function is changed.

## Value

`f`, with the three bindings prepended to its body when it is a function
that needs them, and unchanged otherwise.

## Details

RTMB's replacements are LEXICAL. frmtmb's own objective code binds them
where it needs them (`"c" <- RTMB::ADoverload("c")`), but a function the
USER wrote inherits nothing from those bindings, and base
[`c()`](https://rdrr.io/r/base/c.html) silently strips the advector
class when a plain numeric comes first. The result is a wrong gradient
with no error.

A body rewrite, not an environment injection, on purpose: the byte-code
compiler may inline base builtins unless the shadowing is visible in the
body it compiles, so an environment binding is not reliably honored
while these assignments always are. `ADoverload()` resolves at CALL
time - base versions off the tape, AD versions on it - so the wrapped
function stays correct in both worlds, and the cost is three lookups per
call of a function that is taped once.

A function that already names `ADoverload` manages its own scope and is
returned untouched, so a hand-placed binding is respected and a double
wrap never stacks. Helpers the wrapped function CALLS still use their
own environments - lexical scope cannot be injected transitively - so a
helper that builds advectors needs the explicit spelling of its own.

## See also

[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
and
[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md),
whose function slots frmtmb wraps with this

## Examples

``` r
dyn <- function(t, y, p) list(c(-p[1] * y[1], p[1] * y[1]))
body(frmtmb_ad_overload(dyn))[[2]]
#> "c" <- RTMB::ADoverload("c")
```
