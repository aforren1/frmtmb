# Reviewer 2, item 5: compare the objective and gradient, lane against
# lib7 (the all-rows build) and against base on the plain design.
L <- readRDS("dev/phase3b-review2/rows-lane.rds")
S <- readRDS("dev/phase3b-review2/rows-lib7.rds")
B <- readRDS("dev/phase3b-review2/rows-base.rds")
for (nm in names(L)) {
  rf <- max(abs(L[[nm]]$fn - S[[nm]]$fn) / abs(S[[nm]]$fn))
  rg <- max(abs(L[[nm]]$gr - S[[nm]]$gr) / pmax(abs(S[[nm]]$gr), 1))
  cat(sprintf("%-12s lane vs lib7: fn rel %.3e identical %s; grad %.3e  codes %s\n", nm, rf,
              identical(L[[nm]]$fn, S[[nm]]$fn), rg,
              paste(names(L[[nm]]$codes), L[[nm]]$codes, sep = ":", collapse = " ")))
}
cat(sprintf("none_cens0 (lane) vs none_plain (base): fn identical %s, rel %.3e\n",
            identical(L$none_cens0$fn, B$none_plain$fn),
            max(abs(L$none_cens0$fn - B$none_plain$fn) / abs(B$none_plain$fn))))
cat(sprintf("none_plain lane vs base: identical %s\n", identical(L$none_plain$fn, B$none_plain$fn)))
