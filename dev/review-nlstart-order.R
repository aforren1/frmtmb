# Review of lane nlstart: the same model, the same data, its two terms
# written in the other order (coefficients renamed to match), from the
# default start. The probe's 4 patterns depend on the coefficient order.
# Usage: REV_ARM=after|before Rscript dev/review-nlstart-order.R
source("dev/review-nlstart-common.R")
src <- readLines("dev/review-nlstart-fresh.R")
eval(parse(text = src[grep("^F <- list", src):(grep("^cat\\(\"ROWHEAD", src) - 1L)]))
rows <- list(c("bernoulli", 33), c("bernoulli", 22), c("binomial", 33),
             c("poisson", 11), c("gamma", 33), c("cumulative", 22))
for (r in rows) {
  fam <- r[1]; seed <- as.numeric(r[2])
  d <- gen(fam, "two_sat_rev", seed * 1000 + match("two_sat_rev", names(B)))
  lhs <- if (fam == "binomial") "y | trials(nt)" else "y"
  f1 <- bf(as.formula(paste(lhs, "~ b1 * (1 - exp(-b2 * z)) + b3 * (1 - exp(-b4 * x))")),
           b1 + b2 + b3 + b4 ~ 1, nl = TRUE)
  f2 <- bf(as.formula(paste(lhs, "~ b3 * (1 - exp(-b4 * x)) + b1 * (1 - exp(-b2 * z))")),
           b3 + b4 + b1 + b2 ~ 1, nl = TRUE)
  f3 <- bf(as.formula(paste(lhs, "~ c1 * (1 - exp(-c2 * x)) + c3 * (1 - exp(-c4 * z))")),
           c1 + c2 + c3 + c4 ~ 1, nl = TRUE)
  a1 <- tt(frm(f1, family = F[[fam]], data = d))
  a2 <- tt(frm(f2, family = F[[fam]], data = d))
  a3 <- tt(frm(f3, family = F[[fam]], data = d))
  cat(sprintf("%-10s %d  b1..b4 (z term first): %.4f w%d | same, nlpars declared b3+b4+b1+b2: %.4f w%d | renamed c1..c4 (x term first): %.4f w%d\n",
              fam, seed, llv(a1), length(attr(a1, "w")), llv(a2),
              length(attr(a2, "w")), llv(a3), length(attr(a3, "w"))))
}
