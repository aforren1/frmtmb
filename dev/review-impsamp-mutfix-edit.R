# Second edit of the fix_h copy (outside the worktree): the joint Hessian's
# diagonal, read once before the loop of imp_anchor_slopes(). Written with the
# Write tool because bash heredocs eat backslashes here.
p <- "C:/Users/adf44/source/r/rev-impsamp-src/fix_h/R/importance.R"
s <- readLines(p)
i <- grep("for (j in seq_along(fix)) {", s, fixed = TRUE)
stopifnot(length(i) == 1L)
s <- append(s, "  djoint <- as.numeric(Matrix::diag(env$spHess(full, random = FALSE)))[fix]",
            i - 1L)
writeLines(s, p)
cat("edited at line", i, "\n")
