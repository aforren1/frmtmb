# Extracted from test-parse.R:58

# test -------------------------------------------------------------------------
spec <- frm(bf(y | trials(n) ~ x) + binomial(),
                 data = NULL, dry_run = "spec")
expect_named(spec$responses[[1]]$aterms, "trials")
expect_error(bf_spec <- frm(bf(y | se(s2) ~ x) + gaussian(),
                                 data = NULL, dry_run = "spec"),
               "se")
