# Find `$` reads of a partial-matching container

A testing aid for frmtmb and its extensions, not part of a fitting
workflow. It reports every place in a package where `$` is used on one
of the containers whose slot names collide under partial matching, so an
extension can assert in its own test suite that it has none.

## Usage

``` r
frm_hazard_reads(package)
```

## Arguments

- package:

  A package name, or an environment holding functions (a namespace,
  which is what a package name resolves to).

## Value

A character vector, one entry per read, spelled
`"function: container$slot"`, sorted by function name. A package with no
hazard read gives `character(0)`.

## Details

`$` on a list partial-matches: when no name matches exactly, R returns
the slot whose name the requested one is a prefix of, rather than
`NULL`. frmtmb keeps several containers whose vocabulary is open (a
distributional parameter is named by the user) or prefix-colliding by
construction (`theta` against `thetaac`, `b` against `beta`, `se`
against `se_sigma`), and a `$` read of an absent slot on one of those
silently returns a neighbor. `[[` returns `NULL` instead, which is what
such call sites assume they are getting.

The rule is keyed on the container's NAME, not on `$` in general:
`fit$obj`, `x$fit` and `object$draws` are untouched, and a chain is
judged one link at a time by its own left-hand name, which is why
`fit$frame[["par_template"]]` is right and `fit$frame$par_template` is
not. The names are therefore RESERVED: binding one of them to something
that is not the container it names is itself a hit, and the fix is to
rename the local rather than to exempt it.

The scan reads the abstract syntax tree of every function in the
namespace, so it sees `$` in code and not in a string or a comment, and
it works against an installed package, where the sources are no longer
on disk. What it cannot see is top-level package code that is not inside
a function, because only that code's result is installed.

## Examples

``` r
# in a package's own test suite:
#   expect_identical(frmtmb::frm_hazard_reads("frmtmb.spline"),
#                    character(0))
frm_hazard_reads("frmtmb")
#> character(0)
```
