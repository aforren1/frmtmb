# Lane ordstart: what the "big but not flagged" fits are. cum_sratio
# seed 16, default fit: which outer parameter is large, and how the
# log-likelihood moves when that parameter is pushed further out.
#   Rscript dev/ordstart-probe16.R [family] [seed]
args <- commandArgs(trailingOnly = TRUE)
famname <- if (length(args)) args[1] else "cum_sratio"
s <- if (length(args) > 1) as.integer(args[2]) else 16L
.libPaths(c("C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
fam <- switch(famname, cum_cum = mixture(cumulative(), cumulative()),
              cum_sratio = mixture(cumulative(), sratio()))
d <- mk(s, mix = TRUE)
print(table(d$y))
f <- frm(bf(y ~ x + z), family = fam, data = d)
p <- f$opt$par
print(round(p, 3))
print(frmtmb:::mixture_ord_degeneracy(f, names(f$spec$responses)[1]))
i <- which.max(abs(p))
cat("largest:", names(p)[i], "index", i, "value", p[i], "\n")
for (dlt in c(-10, -5, -1, 0, 1, 5, 10)) {
  q <- p
  q[i] <- q[i] + dlt * sign(p[i])
  cat(sprintf("shift %+4g: nll %.10f  grad %.3g\n", dlt, f$obj$fn(q),
              f$obj$gr(q)[i]))
}
print(f$obj$gr(p))
