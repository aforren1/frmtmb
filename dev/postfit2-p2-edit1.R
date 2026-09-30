# One-off source edit (punch round 2): remove the four dead new-level
# functions of the round-0 construction from R/conditional-effects.R.
f <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2/R/conditional-effects.R"
s <- readLines(f)
a <- grep("^#' WHICH GROUP a bootstrap's draws belong to, as a comparable value[.]$", s)
b <- grep("^ce_draw_new_levels <- function", s)
stopifnot(length(a) == 1, length(b) == 1, b > a)
e <- b + which(s[(b + 1):length(s)] == "}")[1]
stopifnot(s[e + 1] == "")
s <- s[-(a:(e + 1))]
writeLines(s, f)
cat("removed", e + 1 - a + 1, "lines\n")
