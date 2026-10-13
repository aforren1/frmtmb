# Lane optrap: test-se-check.R's one-level group with a bound on x
# (`y ~ x + (1 | g1)`, g1 one level, ub = 0.1 on b), over seeds, with
# and without the bound: how often each build ends at a nonzero code.
#   Rscript dev/optrap-onelevel.R <lib or ""> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 101:140
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
X <- list()
for (s in seeds) for (pr in c(TRUE, FALSE)) {
  set.seed(s)
  d3 <- data.frame(x = rnorm(300), g1 = factor("only"))
  d3$y <- rnorm(300, 1 + 2 * d3$x, 1)
  f <- suppressWarnings(frm(bf(y ~ x + (1 | g1)), family = gaussian(),
    data = d3, prior = if (pr) set_prior("", class = "b", ub = 0.1)))
  X[[length(X) + 1L]] <- data.frame(seed = s, bound = pr,
                                    code = f$opt$convergence,
                                    obj = f$opt$objective,
                                    evals = f$opt$evaluations[1])
}
X <- do.call(rbind, X)
write.table(X, sprintf("dev/optrap-log/onelevel-%s.tsv",
                       basename(dirname(find.package("frmtmb")))),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(table(bound = X$bound, code = X$code))
