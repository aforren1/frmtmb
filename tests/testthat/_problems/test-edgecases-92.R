# Extracted from test-edgecases.R:92

# test -------------------------------------------------------------------------
dd <- data.frame(y = c(2, 3, 5), n = c(5, 5, 4), x = 1:3)
expect_error(frm(bf(y | trials(n) ~ x) + binomial(), data = dd),
               "\\[0, trials\\]")
dd2 <- data.frame(y = c(0.5, 1), n = c(2, 2), x = 1:2)
expect_error(frm(bf(y | trials(n) ~ x) + binomial(), data = dd2),
               "integer")
