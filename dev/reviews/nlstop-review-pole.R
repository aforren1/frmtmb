# Review of lane nlstop: does any default or user start of Vm * xp / (K +
# xp) end at a pole (K in [-max(xp), -min(xp)]) with no warning? Rows of
# the review's fresh set (dev/nlstop-fresh.R generator): bernoulli 55 and
# binomial 11, 44. Starts: a grid of (Vm, K) around and between the poles.
# Usage: Rscript dev/reviews/nlstop-review-pole.R
source("dev/nlstop-common.R")
src <- parse("dev/nlstop-fresh.R")
keep <- vapply(src, function(e) {
  is.call(e) && identical(e[[1L]], as.name("<-")) &&
    as.character(e[[2L]])[1] %in% c("F", "pos", "B", "gen", "form")
}, NA)
for (e in src[keep]) eval(e)
cat("ROWHEAD fam seed start ll code wclass Vm K pole\n")
for (cs in list(c("bernoulli", 55), c("binomial", 11), c("binomial", 44))) {
  fam <- cs[1]; seed <- as.numeric(cs[2])
  d <- gen(fam, "mm", seed * 1000 + match("mm", names(B)))
  fm <- form(fam, "mm")
  rx <- range(d$xp)
  starts <- list(default = NULL)
  for (K in c(-6, -5.05, -4.95, -2.5, -0.15, -0.05, 0.5, 50)) {
    for (Vm in c(-2, 0.5, 5)) {
      starts[[sprintf("Vm%g_K%g", Vm, K)]] <-
        list(beta = c(Vm, K))
    }
  }
  for (sn in names(starts)) {
    r <- tt(frm(fm, family = F[[fam]], data = d, start = starts[[sn]],
                se = TRUE))
    if (iserr(r)) {
      cat(sprintf("ROW %s %g %s ERR %s\n", fam, seed, sn,
                  gsub("\\s+", " ", substr(r$msg, 1, 100))))
      next
    }
    b <- unname(r$estimates$beta)
    pole <- b[2] <= -rx[1] + 1e-8 && b[2] >= -rx[2] - 1e-8
    cat(sprintf("ROW %s %g %s %.6f %s %s %.5g %.6g %s\n", fam, seed, sn,
                llv(r), r$opt$convergence, wsum(r), b[1], b[2],
                if (pole) "POLE" else "-"))
    flush.console()
  }
}
cat("DONE\n")
