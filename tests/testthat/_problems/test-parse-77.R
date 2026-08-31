# Extracted from test-parse.R:77

# test -------------------------------------------------------------------------
expect_error(frm(bf(y ~ rr(x | g)) + gaussian(),
                      data = NULL, dry_run = "spec"),
               "not supported yet")
