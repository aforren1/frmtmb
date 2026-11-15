# Review of lane silent2: 5353 contrasts, the cases silent2-review-core.R
# did not reach. Usage: Rscript dev/reviews/silent2-review-contr.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-silent2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", find.package("frmtmb"), "\n")
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
set.seed(7)
dc <- data.frame(f = factor(sample(c("a", "b", "c"), 90, TRUE)),
                 g = factor(rep(1:10, each = 9)), x = rnorm(90),
                 cyl = sample(c(4, 6, 8), 90, TRUE))
dc$byv <- c("p", "q")[1 + (as.integer(dc$g) %% 2)]
dc$g2 <- factor(rep(c(2:10, 1), each = 9), levels = 1:10)
dc$y <- as.numeric(dc$f) + rnorm(10)[dc$g] + rnorm(90)
nd <- dc[c(1:6, 30), ]
cmp <- function(lab, fml, set_opt, reset_opt, ...) {
  r1 <- withr::with_options(list(contrasts = set_opt), {
    f <- q(fml, data = dc, ...)
    list(f = f, p = tryCatch(fitted(f, newdata = nd)[, "Estimate"],
                             error = function(e) conditionMessage(e)),
         ins = fitted(f)[c(1:6, 30), "Estimate"])
  })
  p2 <- withr::with_options(list(contrasts = reset_opt),
          tryCatch(fitted(r1$f, newdata = nd)[, "Estimate"],
                   error = function(e) conditionMessage(e)))
  if (is.character(r1$p) || is.character(p2)) {
    cat(sprintf("%-58s ERROR fit-opt: %s | reset: %s\n", lab,
                substr(if (is.character(r1$p)) r1$p else "ok", 1, 90),
                substr(if (is.character(p2)) p2 else "ok", 1, 90)))
  } else {
    cat(sprintf("%-58s fit-opt vs in-sample %.3g ; reset vs fit-opt %.3g\n",
                lab, max(abs(unname(r1$p) - unname(r1$ins))),
                max(abs(unname(p2) - unname(r1$p)))))
  }
}
sumc <- c("contr.sum", "contr.poly")
trc <- c("contr.treatment", "contr.poly")
helc <- c("contr.helmert", "contr.poly")
cmp("(1+f|g) sum -> helmert (same column names)", y ~ f + (1 + f | g), sumc, helc)
cmp("(1+factor(cyl)|g) sum -> treatment", y ~ factor(cyl) + (1 + factor(cyl) | g), sumc, trc)
cmp("(1+factor(cyl)|g) sum -> helmert (same names)", y ~ factor(cyl) + (1 + factor(cyl) | g), sumc, helc)
cmp("(1+f|gr(g, by=byv)) sum -> treatment", y ~ f + (1 + f | gr(g, by = byv)), sumc, trc)
cmp("sigma ~ (1+f|g) sum -> treatment", bf(y ~ f, sigma ~ (1 + f | g)), sumc, trc)
cmp("(1+f|g) helmert -> treatment", y ~ f + (1 + f | g), helc, trc)
cmp("(1+f|mm(g,g2)) sum -> helmert", y ~ x + (1 + f | mm(g, g2)), sumc, helc)
cmp("(1+f|g) + (1|g2) partial re_formula handled elsewhere", y ~ f + (1 + f | g) + (1 | g2), sumc, trc)
# factor's own contrasts attribute at the fit, newdata without it
dca <- dc
contrasts(dca$f) <- contr.sum(3)
fa <- q(y ~ f + (1 + f | g), data = dca)
nd2 <- nd
nd2$f <- factor(as.character(nd2$f), levels = levels(dc$f))
pa <- tryCatch(fitted(fa, newdata = nd2)[, "Estimate"], error = function(e) conditionMessage(e))
cat(sprintf("%-58s %s\n", "data attribute contr.sum, newdata plain",
            if (is.character(pa)) paste("ERROR", substr(pa, 1, 120)) else
              format(max(abs(unname(pa) - unname(fitted(fa)[c(1:6, 30), "Estimate"]))), digits = 3)))
# re_formula naming part of a term, after reset
withr::with_options(list(contrasts = sumc), {
  fp <- q(y ~ f + (1 + f | g) + (1 | g2), data = dc)
  r0 <- fitted(fp, newdata = nd, re_formula = ~ (1 + f | g))[, "Estimate"]
  r0b <- tryCatch(fitted(fp, newdata = nd, re_formula = ~ (1 | g))[, "Estimate"],
                  error = function(e) conditionMessage(e))
})
r1 <- tryCatch(fitted(fp, newdata = nd, re_formula = ~ (1 + f | g))[, "Estimate"],
               error = function(e) conditionMessage(e))
cat(sprintf("%-58s %s\n", "re_formula ~(1+f|g) after reset",
            if (is.character(r1)) paste("ERROR", substr(r1, 1, 120)) else
              format(max(abs(r1 - r0)), digits = 3)))
r1b <- tryCatch(fitted(fp, newdata = nd, re_formula = ~ (1 | g))[, "Estimate"],
                error = function(e) conditionMessage(e))
cat(sprintf("%-58s %s\n", "re_formula ~(1|g) (part) after reset",
            if (is.character(r1b) || is.character(r0b)) paste("ERROR/at-fit:", substr(paste(r1b), 1, 80), "|", substr(paste(r0b), 1, 60)) else
              format(max(abs(r1b - r0b)), digits = 3)))
cat("DONE\n")
