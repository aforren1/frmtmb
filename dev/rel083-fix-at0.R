# Merge step 2 (famend x priors2), test-prior-table-brms.R's one conflict:
# resolved with priors2's corrected note (its findings, fix round 2,
# "Merge note correction"), which its short check's scratch merge ran
# (105 pass): a caller's start at the box's end reaches m8's mode, code 0.
# famend's code-1 caller-start block and priors2's cold-start at0 (which
# famend's loop over sp() and sp(lb = 0, ub = 1) already covers) go.
f <- "tests/testthat/test-prior-table-brms.R"
x <- readLines(f)
crlf <- any(grepl("\r$", readLines(f, warn = FALSE)))
a <- grep("^<<<<<<< rel", x); m <- grep("^=======", x); b <- grep("^>>>>>>> famend", x)
stopifnot(length(a) == 1, length(b) == 1, length(m) == 1, a < m, m < b)
new <- c(
"  # the caller's start at the box's end: 0.82.0 stopped on it at an",
"  # infinite objective (code 1); the open chart of the end",
"  # (open_end_coords()) reaches the mode lb = 1e-8 finds",
"  at0 <- frm(yneg ~ z, data = d, prior = sp(lb = 0, ub = 1),",
"             start = list(beta = c(1, 0)))",
"  expect_identical(at0$opt$convergence, 0L)",
"  expect_lte(abs(at0$opt$objective - m8$opt$objective),",
"             1e-8 * abs(m8$opt$objective))")
x <- c(x[seq_len(a - 1)], new, x[-seq_len(b)])
con <- file(f, "wb"); writeLines(x, con, sep = "\r\n"); close(con)
cat("resolved\n")
