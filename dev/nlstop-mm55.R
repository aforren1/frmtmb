# Lane nlstop, items 4960 and 4902: bernoulli Michaelis-Menten seed 55
# of the review's fresh set. The returned fit, its standard errors, the
# objective along the ridge (Vm, K) -> s (Vm, K), and the point the
# reference start reached (K = -5) against the data's xp range.
# Usage: NLSTOP_ARM=before|after Rscript dev/nlstop-mm55.R [fam] [seed]
a <- commandArgs(TRUE)
fam <- if (length(a)) a[1] else "bernoulli"
seed <- if (length(a) > 1) as.numeric(a[2]) else 55
source("dev/nlstop-common.R")
src <- parse("dev/nlstop-fresh.R")
keep <- vapply(src, function(e) {
  is.call(e) && identical(e[[1L]], as.name("<-")) &&
    as.character(e[[2L]])[1] %in% c("F", "pos", "B", "gen", "form")
}, NA)
for (e in src[keep]) eval(e)
d <- gen(fam, "mm", seed * 1000 + match("mm", names(B)))
fm <- form(fam, "mm")
cat("xp range:", range(d$xp), " rows at max xp:", sum(d$xp == max(d$xp)),
    " y there:", d$y[d$xp == max(d$xp)], "\n")
r <- tt(frm(fm, family = F[[fam]], data = d))
cat("default ll", format(llv(r), digits = 10), "code", r$opt$convergence,
    "warnings", length(attr(r, "w")), "\n")
print(suppressWarnings(fixef(r)))
p <- r$opt$par
cat("par:", format(p, digits = 8), "\n")
for (s in c(0.5, 2, 10, 100, 1e4)) {
  q <- p; q[1:2] <- q[1:2] * s
  cat(sprintf("ridge s = %g: ll change %.3e\n", s,
              -(r$obj$fn(q) - r$opt$objective)))
}
lin <- tt(frm(bf(y ~ 0 + xp), family = F[[fam]], data = d))
cat("linear limit y ~ 0 + xp: ll", format(llv(lin), digits = 10),
    "slope", fixef(lin)[1, 1], " Vm/K at the fit", p[1] / p[2], "\n")
