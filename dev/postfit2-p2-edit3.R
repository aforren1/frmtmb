# One-off source edit: the re_formula reason is checked after the
# simulation one, which says more when both differ.
f <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2/R/conditional-effects.R"
s <- readLines(f)
a <- grep("} else if (!is.null(boot$ce_re) &&", s, fixed = TRUE)
b <- grep("} else if (!identical(boot$ce_kept %||% integer(0), kept)) {", s, fixed = TRUE)
stopifnot(length(a) == 1, b == a + 6)
re <- s[a:(a + 5)]
e <- b + which(grepl("^               } else \{$", s[(b + 1):length(s)]))[1]
kp <- s[b:(e - 1)]
s <- c(s[seq_len(a - 1)], kp, re, s[e:length(s)])
writeLines(s, f)
cat("ok\n")
