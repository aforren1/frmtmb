# Lane sampfix nits round (S3, M1): the frm_sample() doc paragraphs.
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/extensions/frmtmb.sample/R/sample.R"
s <- gsub("\r", "", readLines(f))
a <- which(s == "#' *Why they refuse rather than fill the values in.* Each draw's")
b <- which(s == "#' to get these quantities.")
stopifnot(length(a) == 1L, length(b) == 1L, b > a)
new <- c(
"#' *Why they refuse rather than fill the values in.* `laplace = TRUE`",
"#' samples an approximate posterior, and its use is to check that",
"#' approximation against the full one (see [check_laplace()]). For",
"#' predictions, and for anything else that conditions on the group",
"#' effects, sample without `laplace = TRUE`. Filling the values in would",
"#' also be wrong or random: a prediction at each draw's conditional",
"#' modes leaves out the uncertainty of the group effects given the",
"#' outer parameters, and understated the spread of the in-sample",
"#' `posterior_epred()` by up to a factor of four on a gaussian model",
"#' with 12 groups of 10 rows; drawing them from their Laplace",
"#' conditional matched the full draws there, but it would make",
"#' `posterior_epred()` depend on the random seed and disagree with",
"#' `log_lik()` and `ranef()` about what a draw contains.")
s <- c(s[seq_len(a - 1L)], new, s[(b + 1L):length(s)])
old <- c(
"#' refusal names the function and says what to do instead. Whether a",
"#' call reads an integrated value is found out by evaluating it at one",
"#' draw with those values set to `NA` and to `0`: a result that changes",
"#' read them.")
i <- which(s == old[1L])
stopifnot(length(i) == 1L, identical(s[i + 0:3], old))
s <- c(s[seq_len(i - 1L)],
"#' refusal names the function you called and says what to do instead;",
"#' it suggests `re_formula = NA` only where that function takes it and",
"#' it computes on the model. Whether a call reads an integrated value is",
"#' found out by evaluating it at one draw with those values set to `NA`",
"#' and to `0`: a result that changes read them. A draw whose result has",
"#' non-finite values where the first draw's has none is refused the same",
"#' way, so a value that vanishes at the first draw only cannot slip",
"#' through as NaN.",
       s[(i + 4L):length(s)])
writeLines(s, f)
cat("ok\n")
