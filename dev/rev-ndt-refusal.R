# rev-ndt: does the new "bound is not set yet" refusal FALSE ALARM?
#
# A check that fires on a correct model is worse than no check, so the
# rate is measured on models a user would actually write rather than
# asserted. Every case below is labelled with what SHOULD happen, and
# both builds are run, because a refusal that 0.6.0 also gave is not a
# cost this change introduced.
#
# REV_ARM = "new" | "old". Seed 51.

arm <- Sys.getenv("REV_ARM", "new")
.libPaths(if (identical(arm, "new")) {
  c("C:/Users/adf44/source/r/rev-ndt-lib",
    "C:/Users/adf44/source/r/reflib-r2",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
} else {
  c("C:/Users/adf44/source/r/reflib-r2",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
})
suppressMessages({library(frmtmb); library(frmtmb.eam)})

set.seed(51)
d <- ddm_simulate(250, mu = 1.0, bs = 1.4, ndt = 0.30, bias = 0.5)
k <- sample(250, 15)
d$rt[k] <- runif(15, 0.12, 2.5)
d$g <- factor(rep(c("a", "b"), length.out = nrow(d)))
lo <- min(d$rt)

say <- function(label, want, expr) {
  r <- tryCatch({force(expr); "OK"},
                error = function(e) paste0("REFUSED: ",
                  substr(conditionMessage(e), 1, 90)))
  hit <- if (identical(r, "OK")) "fit" else "refused"
  flag <- if (identical(hit, want)) "  " else "**"
  cat(sprintf("%s %-46s want=%-7s got=%-7s %s\n", flag, label, want,
              hit, if (identical(r, "OK")) "" else r))
  identical(hit, want)
}

res <- logical(0)
mx <- function(...) frmtmb::mixture(...)

cat("arm", arm, " min(rt)", sprintf("%.4f", lo), "\n")

# ---- models that MUST fit -------------------------------------------
res <- c(res, say("plain wiener(), no mixture", "fit",
  frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
      family = wiener(), data = d)))
res <- c(res, say("wiener() with a grouping in the data", "fit",
  frm(bf(rt | dec(upper) ~ g, bs ~ 1, ndt ~ 1, bias = 0.5),
      family = wiener(), data = d)))
res <- c(res, say("mixture(wiener(max_ndt below min rt), lognormal)",
  "fit",
  frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
      family = mx(wiener(max_ndt = lo * 0.9), frmtmb::lognormal()),
      data = d)))
res <- c(res, say("mixture(wiener(max_ndt, allow_unreachable), lnorm)",
  "fit",
  frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
      family = mx(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
                  frmtmb::lognormal()), data = d)))
res <- c(res, say("mixture of TWO wiener(max_ndt)", "fit",
  frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5, bias2 = 0.5),
      family = mx(wiener(max_ndt = lo * 0.9),
                  wiener(max_ndt = lo * 0.5)), data = d)))
# frm_simulate() runs assemble_frame(), so the family IS finalized on
# this path and the bound is set from the response. The parameter shape
# is read off a dry run so that the case tests the family and not the
# caller's arithmetic.
res <- c(res, say("frm_simulate() from a bare wiener()", "fit", {
  fsim <- bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5)
  set.seed(303)
  sm <- frm_simulate(fsim, data = d, family = wiener(),
                     newparams = list(beta = 1,
                                      betad = c(log(1.4), 0, 0)),
                     nsim = 1)
  cat("   frm_simulate draw checksum",
      sprintf("%.17g", sum(as.numeric(sm[["rt"]]))), "
")
  sm
}))

# rdm() and lba() cannot go inside mixture() at all: frmtmb's
# mixture() requires every component to have a `mu` dpar and neither
# race family declares one. The case is unwritable rather than
# refused, on both builds.

# ---- models that MUST be refused ------------------------------------
res <- c(res, say("mixture(bare wiener(), lognormal)", "refused",
  frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
      family = mx(wiener(), frmtmb::lognormal()), data = d)))
res <- c(res, say("ndt_group() inside a mixture", "refused",
  frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bias1 = 0.5),
      family = mx(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
                  frmtmb::lognormal()), data = d)))
res <- c(res, say("max_ndt and ndt_group() together", "refused",
  frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
         bias = 0.5), family = wiener(max_ndt = lo * 0.9), data = d)))

cat(sprintf("\narm %s: %d of %d cases behaved as wanted\n", arm,
            sum(res), length(res)))
