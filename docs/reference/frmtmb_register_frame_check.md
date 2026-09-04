# Check an assembled model frame from another package

Registers a function that frmtmb calls on every model frame it
assembles, after the response and the design matrices are built and
before the objective is taped. A package that adds model syntax of its
own uses this to refuse a data problem that only the assembled frame can
show. Register from the contributing package's `.onLoad()`.

## Usage

``` r
frmtmb_register_frame_check(fn)
```

## Arguments

- fn:

  A function of two arguments, the model specification and the assembled
  frame.

## Value

`NULL`, invisibly. Called for the registration.

## Details

The check is called as `fn(spec, frame)` and its return value is
discarded. To refuse a model, call
[`stop()`](https://rdrr.io/r/base/stop.html) with a message that says
what is wrong and what the user can do instead. A check with nothing to
say returns without error. The frame is a read-only view: a check must
not try to change the model.

`frame` is a named list. These elements are the contract:

- `spec`:

  The parsed model specification, the same object the first argument
  carries.

- `n_obs`:

  Number of observations.

- `y`:

  The response, already coerced.

- `y_levels`:

  Response levels for a categorical family, else `NULL`.

- `aterm_values`:

  Values of the addition terms, by response.

- `data_frame`:

  The model frame, one row per observation.

- `linpred`:

  `function(resp_name, dpar)` giving the linear predictor of one
  distributional or nonlinear parameter, or `NULL` when there is none.
  Elements `X` and `Z` are its fixed and random design matrices. Use
  this accessor. The key format of the `linpreds` list itself is
  internal and can change.

Other elements of `frame` are internal assembly structures. They are
present, but they are not part of the contract.

## See also

[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
for the family-side protocol,
[`frmtmb_register_aterm()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_aterm.md)
and
[`frmtmb_register_compat()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_compat.md)
for the other two registries, and
[frmtmb-extension-api](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
for the accessors an extension may use

## Examples

``` r
# the shape of a check: refuse what only the assembled frame shows
check_not_constant <- function(spec, frame) {
  if (length(unique(frame$y)) == 1L) {
    stop("The response takes one value, so nothing can be estimated.",
         call. = FALSE)
  }
}

# a registration lasts for the session, so it belongs in the
# contributing package's .onLoad(), not in a script
if (FALSE) { # \dontrun{
frmtmb_register_frame_check(check_not_constant)
} # }
```
