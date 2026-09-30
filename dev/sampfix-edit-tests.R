# Lane sampfix nits round: insert the new laplace tests.
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/extensions/frmtmb.sample/tests/testthat/test-laplace-draws.R"
s <- gsub("\r", "", readLines(f))
i <- which(s == "## ---- frm_sample() itself --------------------------------------------")
stopifnot(length(i) == 1L)
add <- readLines("C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-newtests.txt")
writeLines(c(s[seq_len(i - 1L)], add, s[i:length(s)]), f)
cat("ok\n")
