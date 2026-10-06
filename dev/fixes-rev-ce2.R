# Reviewer of lane fixes, claim 3: a factor read by an addition term,
# against brms 2.23.0 at frmtmb's estimates.
#   Rscript dev/fixes-rev-ce2.R <lib>
src <- readLines("dev/fixes-rev-ce.R")
src <- src[seq_len(grep("^run[(]\"trunc[(]lb = min", src) - 1L)]
eval(parse(text = src))
d$lb2 <- ifelse(d$g == "a", -2.5, -2)
d <- d[d$y > d$lb2, ]
run("trunc(lb = ifelse(g == 'a', -2.5, -2)), g only in the term",
    bf(y | trunc(lb = ifelse(g == "a", -2.5, -2)) ~ x),
    brms::bf(y | trunc(lb = ifelse(g == "a", -2.5, -2)) ~ x))
run("same with conditions g = 'b'",
    bf(y | trunc(lb = ifelse(g == "a", -2.5, -2)) ~ x),
    brms::bf(y | trunc(lb = ifelse(g == "a", -2.5, -2)) ~ x),
    cond = data.frame(g = factor("b", levels = levels(d$g))))
