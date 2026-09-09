# rev-ndt round 2, attack 3: the equivalence claim, at a FIXED
# parameter vector, on all four families that changed plus gddm.
#
# Two fitted optima agreeing is strong but not the claim. The claim is
# that with no ndt_group() nothing moved at all, so the objective and
# its gradient are compared where nothing but arithmetic can differ: at
# the starting vector and at a deliberate off-optimum point.
#
# REV_ARM = "new" | "old". Seed 4242.

arm <- Sys.getenv("REV_ARM", "new")
new_lib <- Sys.getenv("REV_LIB", "C:/Users/adf44/source/r/rev-ndt-lib2")
ref_lib <- "C:/Users/adf44/source/r/reflib-r2"
usr_lib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "new")) c(new_lib, ref_lib, usr_lib)
          else c(ref_lib, usr_lib))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

out <- list()
probe <- function(nm, form, fam, dat) {
  ob <- frm(form, family = fam, data = dat, dry_run = "objective")$obj
  p <- ob$par
  pts <- list(start = p, off = p + seq_along(p) * 0.03)
  for (w in names(pts)) {
    out[[paste0(nm, ".", w, ".fn")]] <<- as.numeric(ob$fn(pts[[w]]))
    out[[paste0(nm, ".", w, ".gr")]] <<- as.numeric(ob$gr(pts[[w]]))
  }
  out[[paste0(nm, ".par")]] <<- as.numeric(p)
  invisible(NULL)
}

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

probe("wiener", bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
      wiener(), d)
probe("wiener_maxndt",
      bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
      wiener(max_ndt = 0.15), d)
probe("wiener_st",
      bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
      wiener(variability = "st"), d)
probe("wiener_svszst",
      bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
      wiener(variability = c("sv", "sz", "st")), d)

set.seed(4242)
dr <- rdm_simulate(300, v = c(2.5, 1.5), A = 0.5, k = 0.5, ndt = 0.2)
probe("rdm", bf(rt | vint(choice) ~ 1, v2 ~ 1, A ~ 1, k ~ 1, ndt ~ 1),
      rdm(2), dr)
dr2 <- dr
dr2$cn <- ifelse(dr2$rt > 0.8, "right", "none")
dr2$rtc <- pmin(dr2$rt, 0.8)
probe("rdm_cens",
      bf(rtc | vint(choice) + cens(cn) ~ 1, v2 ~ 1, A ~ 1, k ~ 1,
         ndt ~ 1), rdm(2), dr2)

set.seed(4242)
dl <- lba_simulate(300, v = c(2.5, 1.5), A = 0.5, k = 0.4, ndt = 0.2)
probe("lba", bf(rt | vint(choice) ~ 1, v2 ~ 1, A ~ 1, k ~ 1, ndt ~ 1),
      lba(2), dl)

set.seed(4242)
dg <- wiener_gng_simulate(400, mu = 1.2, bs = 1.4, ndt = 0.25,
                          deadline = 1.5)
probe("wiener_gng",
      bf(rt | dec(responded) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
      wiener_gng(deadline = 1.5), dg)
probe("wiener_gng_st",
      bf(rt | dec(responded) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
      wiener_gng(deadline = 1.5, variability = "st"), dg)

set.seed(5)
dq <- gddm_simulate(120, mu = 2, bs = 2.5, ndt = 0.25,
                    control = gddm_control(t_max = 2))
dq$cond <- 1L
probe("gddm", bf(rt | vint(upper, cond) ~ 1, bs ~ 1, ndt ~ 1,
                 bias = 0.5),
      gddm(control = gddm_control(t_max = 2, dt = 0.05, ny = 51L)), dq)

saveRDS(out, paste0("C:/Users/adf44/source/r/frmtmb-wt-ndt/dev/",
                    "rev-ndt-density2-", arm, ".rds"))
cat("ARM", arm, "OK", length(out), "probes\n")
