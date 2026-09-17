## Reviewer recheck: positional (level x coefficient) pairing of ranef(ds)
## against the ML ranef() at draw 20; reuses dev/brmsnames-rev2-ranef.R.
##   Rscript dev/brmsnames-rev2-ranef-pos.R
src <- readLines("dev/brmsnames-rev2-ranef.R")
end <- grep("^flat <- function", src) - 1L
eval(parse(text = src[1:end]))
for (nm in c("us", "rr", "ar1", "cs", "toep", "gr_prec", "equalto", "exp", "two_blocks", "smooth_re")) {
  M <- ms[[nm]]
  fit <- q(frm(bf(M[[1]]), family = gaussian(), data = dd, data2 = d2))
  set.seed(4)
  ds <- q(frm_sample(fit, chains = 1, iter = 80, refresh = 0, seed = 4))
  rd <- ranef(ds, summary = FALSE)
  idx <- sns$draws_par_index(fit)
  sh <- sns$draws_fit_at(ds, 20, idx)
  mr <- ranef(sh)
  out <- c()
  for (k in seq_along(rd)) {
    a <- rd[[k]][20, , , drop = FALSE][1, , , drop = TRUE]
    a <- matrix(a, dim(rd[[k]])[2])
    hit <- which(vapply(mr, function(t) all(dim(as.matrix(t)) == dim(a)), TRUE))
    dif <- vapply(hit, function(h) max(abs(as.matrix(mr[[h]]) - a)), 1)
    out <- c(out, sprintf("%s: min positional max|diff| over same-shape ML tables %.3g (ML names %s)",
                          names(rd)[k], min(dif), paste(names(mr)[hit], collapse = ",")))
  }
  cat(sprintf("%-11s %s\n", nm, paste(out, collapse = " | ")))
}
