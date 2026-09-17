# Reviewer recheck round 1: is a collapsed twin identical to the single
# term, bit for bit? Lane build; data as dev/priorform-rev2-special.R.
#   Rscript dev/priorform-rev2-twinll.R
.libPaths(c("C:/Users/adf44/source/r/priorform-lib", "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-priorform/dev/priorform-rev2-special.R")
eval(parse(text = src[grep("^set.seed", src):grep("^d2 <- ", src)]))
pairs <- list(c("y ~ (1 | g/h) + (1 | g/h)", "y ~ (1 | g/h)"),
              c("y ~ (x || g) + (x || g)", "y ~ (x || g)"),
              c("y ~ (1 | p | g) + (1 | p | g)", "y ~ (1 | p | g)"),
              c("y ~ x + (1+x|g) + (1 + x | g)", "y ~ x + (1 + x | g)"),
              c("y ~ x + s(z) * s(z)", "y ~ x + s(z)"),
              c("y ~ (x + s(z))", "y ~ x + s(z)"))
for (p in pairs) {
  f <- lapply(p, function(s) suppressWarnings(suppressMessages(frm(bf(as.formula(s)) + gaussian(), data = d))))
  cat(sprintf("%-32s vs %-22s identical opt: %s  logLik %.10f %.10f\n", p[1], p[2],
              identical(f[[1]]$opt[c("par", "objective")], f[[2]]$opt[c("par", "objective")]),
              as.numeric(logLik(f[[1]])), as.numeric(logLik(f[[2]]))))
}
