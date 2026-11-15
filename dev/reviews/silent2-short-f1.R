# Short check of lane silent2's fix round, F1: transformed factors on the
# left side of a group-level term. Fit under contr.sum, predict under
# contr.sum, contr.helmert and contr.treatment; each prediction against
# the in-sample fitted values of the same rows.
#   Rscript dev/reviews/silent2-short-f1.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-silent2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", find.package("frmtmb"), "\n")
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
set.seed(7)
n <- 120
dc <- data.frame(f = factor(sample(c("a", "b", "c"), n, TRUE)),
                 a = factor(sample(c("p", "q"), n, TRUE)),
                 b = factor(sample(c("u", "v"), n, TRUE)),
                 g = factor(rep(1:12, each = 10)), x = rnorm(n),
                 cyl = sample(c(4, 6, 8), n, TRUE))
dc$gnum <- as.integer(dc$g)
dc$g2num <- c(2:12, 1)[dc$gnum]
dc$byv <- c("p", "q")[1 + (dc$gnum %% 2)]
dc$y <- as.numeric(dc$f) + 0.5 * as.numeric(dc$a) * as.numeric(dc$b) +
  rnorm(12)[dc$g] + 0.3 * dc$x + rnorm(n)
rows <- c(1:6, 30, 55, 90)
sumc <- c("contr.sum", "contr.poly")
opts <- list(sum = sumc, helmert = c("contr.helmert", "contr.poly"),
             treatment = c("contr.treatment", "contr.poly"))
pr <- function(f, nd) {
  tryCatch(unname(fitted(f, newdata = nd)[, "Estimate"]),
           error = function(e) paste("REFUSED:", substr(gsub("\n", " ",
                                     conditionMessage(e)), 1, 110)))
}
cmp <- function(lab, fml, nd_rows = rows, ...) {
  fit <- withr::with_options(list(contrasts = sumc),
                             tryCatch(q(fml, data = dc, ...),
                                      error = function(e) e))
  if (inherits(fit, "error")) {
    cat(sprintf("%-46s FIT ERROR: %s\n", lab,
                substr(conditionMessage(fit), 1, 110)))
    return(invisible())
  }
  ins <- unname(fitted(fit)[nd_rows, "Estimate"])
  nd <- dc[nd_rows, ]
  for (o in names(opts)) {
    p <- withr::with_options(list(contrasts = opts[[o]]), pr(fit, nd))
    res <- if (is.character(p)) p else {
      sprintf("max|diff| %.3g", max(abs(p - ins)))
    }
    cat(sprintf("%-46s %-9s %s\n", lab, o, res))
  }
}
cmp("(1 + factor(cyl):x | g)", y ~ x + (1 + factor(cyl):x | g))
cmp("(1 + interaction(a, b) | g)", y ~ 1 + (1 + interaction(a, b) | g))
cmp("(1 + relevel(f, \"b\") | g)", y ~ 1 + (1 + relevel(f, "b") | g))
cmp("(1 + cut(x, 3) | g), all rows", y ~ 1 + (1 + cut(x, 3) | g),
    nd_rows = seq_len(n))
cmp("(1 + cut(x, 3) | g), 9 rows", y ~ 1 + (1 + cut(x, 3) | g))
cmp("(1 + factor(cyl) | gr(g, by = byv))",
    y ~ 1 + (1 + factor(cyl) | gr(g, by = byv)))
cmp("(1 + f | gr(factor(gnum)))", y ~ 1 + (1 + f | gr(factor(gnum))))
cmp("(1 + factor(cyl) | gr(factor(gnum)))",
    y ~ 1 + (1 + factor(cyl) | gr(factor(gnum))))
cmp("(1 + f | mm(factor(gnum), factor(g2num)))",
    y ~ 1 + (1 + f | mm(factor(gnum), factor(g2num))))
cmp("(1 + factor(cyl) | mm(factor(gnum), factor(g2num)))",
    y ~ 1 + (1 + factor(cyl) | mm(factor(gnum), factor(g2num))))
dc$gm1 <- dc$g
dc$gm2 <- factor(dc$g2num, levels = 1:12)
cmp("(1 + factor(cyl) | mm(gm1, gm2))",
    y ~ 1 + (1 + factor(cyl) | mm(gm1, gm2)))
cmp("(1 + interaction(a, b) | mm(gm1, gm2))",
    y ~ 1 + (1 + interaction(a, b) | mm(gm1, gm2)))
cmp("(1 + factor(cyl) || g)", y ~ 1 + (1 + factor(cyl) || g))
cmp("sigma ~ (1 + factor(cyl) | g)",
    bf(y ~ 1 + (1 | g), sigma ~ 1 + (1 + factor(cyl) | g)))
cat("DONE\n")
