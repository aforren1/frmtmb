# Lane nlstop2, item 5079: debug of nl_pole_escape() on a fit from a
# user start (bernoulli mm 55, start Vm -2, K -2.5).
# Usage: Rscript dev/nlstop2-poledbg.R [fam] [seed] [Vm] [K]
a <- commandArgs(TRUE)
fam <- if (length(a)) a[1] else "bernoulli"
seed <- if (length(a) > 1) as.numeric(a[2]) else 55
st <- if (length(a) > 3) as.numeric(a[3:4]) else c(-2, -2.5)
source("dev/nlstop2-common.R")
fresh_defs()
ns <- asNamespace("frmtmb")
d <- gen(fam, "mm", seed * 1000 + match("mm", names(B)))
r <- tt(frm(form(fam, "mm"), family = F[[fam]], data = d,
            start = list(beta = st)))
cat("fit", llv(r), r$opt$convergence, wsum(r), msum(r),
    format(r$estimates$beta, digits = 6), "\n")
lp <- Filter(function(q) !is.null(q$nl_body), r$frame$linpreds)[[1]]
dens <- ns$nl_denoms(lp$nl_body)
print(dens)
pn <- names(r$obj$par)
for (at in 1:2) {
  wl <- ns$nl_pole_walls(r$frame, lp, r$estimates, "beta", at, dens)
  str(wl)
}
wl <- ns$nl_pole_walls(r$frame, lp, r$estimates, "beta", 2, dens)[[1]]
half <- (wl$hi - wl$lo) / 2
for (v in c(wl$lo - half, wl$hi + half)) {
  p <- r$opt$par
  p[2] <- v
  o <- ns$optimize_obj(r$obj, r$control, ns$fit_outer_box(r), r$par_units,
                       start_par = p)
  cat("from K", v, ": objective", o$objective, "code", o$convergence,
      o$message, "par", format(o$par, digits = 6), "\n")
}
run_used <- function(st, announce) stop("not used")
pe <- ns$nl_pole_escape(function(st, an) {
  frm(form(fam, "mm"), family = F[[fam]], data = d, start = st)
}, r$frame, r)
str(pe[c("ll0", "ll1")])
