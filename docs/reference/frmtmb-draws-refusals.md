# Methods a ported brms script may call that frmtmb does not have

These `brmsfit` methods either describe machinery frmtmb does not use
(Stan code and Stan data) or are brms spellings that have been renamed.
They are defined so that a ported script gets the reason and the
replacement rather than "could not find function", which is what the
vignette-port audit measured most of its post-processing failures as.

## Usage

``` r
stancode(object, ...)

# S3 method for class 'frmtmb_draws'
stancode(object, ...)

standata(object, ...)

# S3 method for class 'frmtmb_draws'
standata(object, ...)

expose_functions(x, ...)

# S3 method for class 'frmtmb_draws'
expose_functions(x, ...)

restructure(x, ...)

# S3 method for class 'frmtmb_draws'
restructure(x, ...)

posterior_samples(x, ...)

# S3 method for class 'frmtmb_draws'
posterior_samples(x, ...)

nsamples(object, ...)

# S3 method for class 'frmtmb_draws'
nsamples(object, ...)

parnames(x, ...)

# S3 method for class 'frmtmb_draws'
parnames(x, ...)
```

## Arguments

- object, x, ...:

  Ignored; these functions always stop.

## Value

These functions never return; they signal an error.
