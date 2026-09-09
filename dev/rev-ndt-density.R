# rev-ndt: is the DENSITY the same function, or only the same algebra?
#
# The backward-compatibility control compares two OPTIMA, which a flat
# ridge can separate for reasons that have nothing to do with the
# change. This compares the objective and its gradient at a FIXED
# parameter vector, where nothing but arithmetic can differ. With no
# ndt_group() the two parameterizations put `ndt` on the same linear
# predictor, so the same `par` means the same model in both arms.
#
# REV_ARM = "new" | "old". Seed 4242.

arm <- Sys.getenv("REV_ARM", "new")
new_lib <- "C:/Users/adf44/source/r/rev-ndt-lib"
ref_lib <- "C:/Users/adf44/source/r/reflib-r2"
usr_lib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "new")) c(new_lib, ref_lib, usr_lib)
          else c(ref_lib, usr_lib))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

set.seed(4242)
ns <- 6L
nt <- 60L
u <- rnorm(ns, 0, 0.12)
s <- rep(seq_len(ns), each = nt)
cond <- rep(rep(0:1, each = nt / 2L), times = ns)
d <- ddm_simulate(ns * nt, mu = 0.4 + 0.9 * cond, bs = 1.4,
                  ndt = 0.25 * exp(u[s]), bias = 0.5)
d$s <- factor(s)
d$cond <- factor(cond, labels = c("a", "b"))

out <- list()
probe <- function(nm, form, fam) {
  ob <- frm(form, family = fam, data = d, dry_run = "objective")$obj
  p <- ob$par
  # a deliberate off-optimum point as well as the start, because the
  # start can sit where the two orderings round the same way
  pts <- list(start = p, off = p + seq_along(p) * 0.03)
  for (w in names(pts)) {
    out[[paste0(nm, ".", w, ".fn")]] <<- as.numeric(ob$fn(pts[[w]]))
    out[[paste0(nm, ".", w, ".gr")]] <<- as.numeric(ob$gr(pts[[w]]))
    out[[paste0(nm, ".", w, ".par")]] <<- as.numeric(pts[[w]])
  }
  invisible(NULL)
}

probe("wiener", bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
      wiener())
probe("wiener_st", bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1,
                      bias = 0.5), wiener(variability = "st"))

set.seed(4242)
dr <- rdm_simulate(300, v = c(2.5, 1.5), A = 0.5, k = 0.5, ndt = 0.2)
obr <- frm(bf(rt | vint(choice) ~ 1, v2 ~ 1, A ~ 1, k ~ 1, ndt ~ 1),
           family = rdm(2), data = dr, dry_run = "objective")$obj
pr <- obr$par
out$rdm.start.fn <- as.numeric(obr$fn(pr))
out$rdm.start.gr <- as.numeric(obr$gr(pr))
out$rdm.off.fn <- as.numeric(obr$fn(pr + seq_along(pr) * 0.03))
out$rdm.off.gr <- as.numeric(obr$gr(pr + seq_along(pr) * 0.03))
out$rdm.start.par <- as.numeric(pr)

saveRDS(out, paste0("C:/Users/adf44/source/r/frmtmb-wt-ndt/dev/",
                    "rev-ndt-density-", arm, ".rds"))
for (k in names(out)) {
  cat(sprintf("%-24s %s\n", k,
              paste(sprintf("%.17g", utils::head(out[[k]], 3)),
                    collapse = " ")))
}
