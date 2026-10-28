## Bug fixes

* `prior =` given a list of `set_prior()` specifications, such as
  `list(sp, sp2)`, is refused with a pointer to `sp + sp2` and
  `c(sp, sp2)`, as brms combines them. It gave a generic message about
  a list of length 2.
