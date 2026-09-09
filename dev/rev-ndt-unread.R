# rev-ndt round 2, attack 6: the frame check that refuses an
# `ndt_group()` no family read.
#
# The guard is `ddm_check_ndt_group_read()`, registered through
# `frmtmb_register_frame_check()`. It refuses when the frame carries an
# `ndt_group` addition-term value and the response's family has no
# per-group floor table, which is the mixture case: `mixture()` unions
# its components' allow-lists and never finalizes them.
#
# Every guard built in the 0.55.1 round failed open on its first try, so
# this runs the case where the guarded thing is ABSENT as well as the
# case where it is present, and then measures the false-alarm rate on
# models a user would actually write.
#
# Seed 51. Worktree build at rev-ndt-lib2.

.libPaths(c(Sys.getenv("REV_LIB",
                       "C:/Users/adf44/source/r/rev-ndt-lib2"),
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

res <- logical(0)
say <- function(lbl, want, expr) {
  r <- tryCatch({force(expr); "fit"},
                error = function(e) paste0("refused|",
                  substr(conditionMessage(e), 1, 78)))
  got <- if (identical(r, "fit")) "fit" else "refused"
  ok <- identical(got, want)
  res <<- c(res, ok)
  cat(sprintf("%s%-46s want=%-7s got=%-7s %s\n", if (ok) "   " else "** ",
              lbl, want, got,
              if (got == "fit") "" else sub("^refused[|]", "", r)))
  invisible(ok)
}

set.seed(51)
d <- ddm_simulate(300, mu = 1.0, bs = 1.4, ndt = 0.30, bias = 0.5)
k <- sample(300, 20)
d$rt[k] <- runif(20, 0.12, 2.5)
d$g <- factor(rep(c("a", "b"), length.out = nrow(d)))
lo <- min(d$rt)
mx <- function(...) frmtmb::mixture(...)

cat("-- the guarded thing PRESENT: ndt_group() nothing reads\n")
say("mixture(wiener(max_ndt), lognormal) + ndt_group", "refused",
    frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bias1 = 0.5),
        family = mx(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
                    frmtmb::lognormal()), data = d))
say("mixture of TWO wiener(max_ndt) + ndt_group", "refused",
    frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bias1 = 0.5,
           bias2 = 0.5),
        family = mx(wiener(max_ndt = lo * 0.9),
                    wiener(max_ndt = lo * 0.5)), data = d))

cat("\n-- the guarded thing ABSENT: every model that DOES read it\n")
say("wiener + ndt_group", "fit",
    frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
           bias = 0.5), family = wiener(), data = d))
say("wiener + ndt_group + weights", "fit", {
  d$w <- 1
  frm(bf(rt | dec(upper) + ndt_group(g) + weights(w) ~ 1, bs ~ 1,
         ndt ~ 1, bias = 0.5), family = wiener(), data = d)
})
say("wiener + ndt_group, random effect on ndt", "fit",
    frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1,
           ndt ~ 1 + (1 | g), bias = 0.5), family = wiener(), data = d))
say("wiener + ndt_group, variability = st", "fit",
    frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
           bias = 0.5), family = wiener(variability = "st"), data = d))
say("rdm + ndt_group", "fit", {
  set.seed(7)
  dr <- rdm_simulate(300, v = c(2.5, 1.5), A = 0.5, k = 0.5, ndt = 0.2)
  dr$g <- factor(rep(c("a", "b"), length.out = nrow(dr)))
  frm(bf(rt | vint(choice) + ndt_group(g) ~ 1, v2 ~ 1, A ~ 1, k ~ 1,
         ndt ~ 1), family = rdm(2), data = dr)
})
say("lba + ndt_group", "fit", {
  set.seed(7)
  dl <- lba_simulate(300, v = c(2.5, 1.5), A = 0.5, k = 0.4, ndt = 0.2)
  dl$g <- factor(rep(c("a", "b"), length.out = nrow(dl)))
  frm(bf(rt | vint(choice) + ndt_group(g) ~ 1, v2 ~ 1, A ~ 1, k ~ 1,
         ndt ~ 1), family = lba(2), data = dl)
})
say("wiener_gng + ndt_group", "fit", {
  set.seed(7)
  dg <- wiener_gng_simulate(400, mu = 1.2, bs = 1.4, ndt = 0.25,
                            deadline = 1.5)
  dg$g <- factor(rep(c("a", "b"), length.out = nrow(dg)))
  frm(bf(rt | dec(responded) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
         bias = 0.5), family = wiener_gng(deadline = 1.5), data = dg)
})
say("frm_simulate() with ndt_group", "fit", {
  f0 <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
               bias = 0.5), family = wiener(), data = d)
  frm_simulate(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
                  bias = 0.5), data = d, family = family(f0),
               newparams = list(beta = 1,
                                betad = c(log(1.4), 0, 0)),
               nsim = 1)
})

cat("\n-- the false-alarm cases: models with no ndt_group() at all\n")
say("plain wiener", "fit",
    frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
        family = wiener(), data = d))
say("mixture(wiener(max_ndt), lognormal), no grouping", "fit",
    frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
        family = mx(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
                    frmtmb::lognormal()), data = d))
say("mixture of two lognormals, no eam family at all", "fit",
    frm(bf(rt ~ 1), family = mx(frmtmb::lognormal(),
                                frmtmb::lognormal()), data = d))
say("a gaussian model with no eam family", "fit",
    frm(rt ~ 1, data = d))
say("gddm, which refuses ndt_group by declaration", "fit", {
  set.seed(5)
  dq <- gddm_simulate(120, mu = 2, bs = 2.5, ndt = 0.25,
                      control = gddm_control(t_max = 2))
  dq$cond <- 1L
  frm(bf(rt | vint(upper, cond) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
      family = gddm(control = gddm_control(t_max = 2, dt = 0.05,
                                           ny = 51L)), data = dq)
})

cat(sprintf("\n%d of %d cases behaved as wanted; false alarms on a ",
            sum(res), length(res)))
cat("correct model: ",
    sum(!res[10:length(res)]), "\n", sep = "")
