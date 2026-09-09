# rev-ndt: edge cases of ndt_time() and of the per-row bound that the
# lane's own tests do not cover.
#
# Seed 4242, worktree build.

.libPaths(c("C:/Users/adf44/source/r/rev-ndt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

say <- function(lbl, expr) {
  r <- tryCatch(expr, error = function(e)
    paste("ERROR:", substr(conditionMessage(e), 1, 100)))
  cat(sprintf("  %-44s %s\n", lbl,
              if (is.character(r)) r
              else paste(sprintf("%.6f", utils::head(as.numeric(r), 3)),
                         collapse = " ")))
  invisible(r)
}

set.seed(4242)
d <- ddm_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.30, bias = 0.5)
d$s <- factor(rep(c("a", "b", "c", "dd"), length.out = nrow(d)))
form <- bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1,
           ndt ~ 1 + (1 | s), bias = 0.5)
fit <- frm(form, family = wiener(), data = d, se = TRUE)
floors <- tapply(d$rt, d$s, min)
cat("group floors:", paste(names(floors),
                           sprintf("%.4f", floors), collapse = " "), "\n")

cat("\n-- ndt_time() basics\n")
say("ndt_time(fit) first 3 rows", ndt_time(fit))
one <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
say("ndt_time(fit, newdata = one per group)", ndt_time(fit, one))
say("ndt_time(fit, re.form = NA)", ndt_time(fit, re.form = NA))
say("ndt_time(fit, newdata=one, re.form = NA)",
    ndt_time(fit, one, re.form = NA))

cat("\n-- the refusals\n")
nd_missing <- one
nd_missing$s <- NULL
say("newdata without the ndt_group column", ndt_time(fit, nd_missing))
nd_new <- one[1L, , drop = FALSE]
nd_new$s <- factor("zz")
say("newdata with a group the fit never saw", ndt_time(fit, nd_new))

cat("\n-- NA rows, where the frame and the prediction differ in length\n")
d2 <- d
d2$x <- rnorm(nrow(d2))
d2$x[c(3L, 17L)] <- NA
f2 <- frm(bf(rt | dec(upper) + ndt_group(s) ~ x, bs ~ 1, ndt ~ 1,
             bias = 0.5), family = wiener(), data = d2)
t2 <- say("ndt_time() with 2 rows dropped by na.action", ndt_time(f2))
if (is.numeric(t2)) {
  cat("     length", length(t2), "of", nrow(d2), " NAs at",
      paste(which(is.na(t2)), collapse = " "), "\n")
}

cat("\n-- a mixture fit, which has no ndt_bound of its own\n")
dm <- d
say("ndt_time() on a mixture fit", {
  fm <- frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
            family = frmtmb::mixture(wiener(max_ndt = 0.3),
                                     frmtmb::lognormal()), data = dm)
  ndt_time(fm)
})

cat("\n-- a gddm fit, where ndt is still a time\n")
set.seed(5)
dq <- gddm_simulate(120, mu = 2, bs = 2.5, ndt = 0.25,
                    control = gddm_control(t_max = 2))
dq$cond <- 1L
fq <- frm(bf(rt | vint(upper, cond) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
          family = gddm(control = gddm_control(t_max = 2, dt = 0.05,
                                               ny = 51L)), data = dq)
say("ndt_time() on a gddm fit", ndt_time(fq))
say("predict(gddm, ndt, response)",
    suppressWarnings(predict(fq, dpar = "ndt", type = "response")))

cat("\n-- does a refit on a SUBSET keep the fitted bound?\n")
# the findings claim the floor table is captured at finalize, so a
# leave-one-out style call on rows that drop a group's fastest trial
# must score against the ORIGINAL bound
k <- which.min(d$rt)
say("ndt_time(fit, newdata = data without its fastest row)",
    ndt_time(fit, d[-k, , drop = FALSE]))
say("  ... the same rows, refitted from scratch", {
  f3 <- frm(form, family = wiener(), data = d[-k, , drop = FALSE])
  ndt_time(f3)
})
