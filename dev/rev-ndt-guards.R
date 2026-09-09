# rev-ndt: the guards, with the guarded thing ABSENT as well as wrong.
#
# Each claim in the compat rows and the help is checked by building the
# case it names. Seed 808, worktree build.

.libPaths(c("C:/Users/adf44/source/r/rev-ndt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

say <- function(lbl, expr) {
  r <- tryCatch({v <- force(expr)
                 paste("OK:", paste(sprintf("%.6f",
                   utils::head(as.numeric(v), 3)), collapse = " "))},
                error = function(e)
                  paste("REFUSED:", substr(conditionMessage(e), 1, 96)))
  cat(sprintf("  %-46s %s\n", lbl, r))
}

# 1. wiener_gng: a group with no GO trials has no bound
cat("\n-- wiener_gng, a group whose trials are all no-go\n")
set.seed(808)
dg <- wiener_gng_simulate(600, mu = 1.2, bs = 1.4, ndt = 0.25,
                          deadline = 1.5)
dg$g <- factor(rep(c("a", "b", "c"), length.out = nrow(dg)))
cat("   responded by group:",
    paste(names(table(dg$g)),
          tapply(dg$responded, dg$g, sum), sep = "=", collapse = " "),
    "\n")
dg2 <- dg
dg2$responded[dg2$g == "c"] <- 0
cat("   after forcing group c to all no-go:",
    paste(names(table(dg2$g)),
          tapply(dg2$responded, dg2$g, sum), sep = "=", collapse = " "),
    "\n")
say("all groups have go trials", {
  f <- frm(bf(rt | dec(responded) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
              bias = 0.5), family = wiener_gng(deadline = 1.5),
           data = dg)
  ndt_time(f)
})
say("group c has none", {
  f <- frm(bf(rt | dec(responded) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
              bias = 0.5), family = wiener_gng(deadline = 1.5),
           data = dg2)
  ndt_time(f)
})
say("and the bound really is over GO rows only", {
  f <- frm(bf(rt | dec(responded) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
              bias = 0.5), family = wiener_gng(deadline = 1.5),
           data = dg)
  bd <- frmtmb::single_response(f)[["family"]][["ndt_bound"]]
  go <- tapply(dg$rt[dg$responded > 0.5], dg$g[dg$responded > 0.5], min)
  all <- tapply(dg$rt, dg$g, min)
  cat("\n     floors    ", paste(sprintf("%.6f", bd$floors),
                                 collapse = " "), "\n")
  cat("     go minima ", paste(sprintf("%.6f", go), collapse = " "),
      "\n")
  cat("     all minima", paste(sprintf("%.6f", all), collapse = " "),
      "\n    ")
  bd$floors
})

# 2. the grouping's type
cat("\n-- how the grouping may be spelled\n")
set.seed(808)
d <- ddm_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
d$gf <- factor(rep(c("a", "b"), length.out = nrow(d)))
d$gc <- as.character(d$gf)
d$gi <- as.integer(d$gf)
d$gl <- d$gi == 1L
d$gn <- as.numeric(d$gi)
d$gna <- d$gi
d$gna[5] <- NA
frm_g <- function(col) {
  form <- as.formula(paste0("rt | dec(upper) + ndt_group(", col, ") ~ 1"))
  frm(bf(form, bs ~ 1, ndt ~ 1, bias = 0.5), family = wiener(),
      data = d)
}
say("factor", ndt_time(frm_g("gf")))
say("integer codes", ndt_time(frm_g("gi")))
say("numeric codes", ndt_time(frm_g("gn")))
say("character (must be refused)", ndt_time(frm_g("gc")))
say("logical (must be refused)", ndt_time(frm_g("gl")))
say("integer with an NA (must be refused)", ndt_time(frm_g("gna")))

# 3. a group of one trial: legal, and the bound is that trial
cat("\n-- a group with a single trial\n")
d2 <- d
d2$g1 <- factor(c("solo", rep("rest", nrow(d) - 1L)))
say("one group has one row", {
  f <- frm(bf(rt | dec(upper) + ndt_group(g1) ~ 1, bs ~ 1, ndt ~ 1,
              bias = 0.5), family = wiener(), data = d2)
  bd <- frmtmb::single_response(f)[["family"]][["ndt_bound"]]
  cat("\n     floors", paste(sprintf("%.6f", bd$floors), collapse = " "),
      " first rt", sprintf("%.6f", d2$rt[1L]), "\n    ")
  ndt_time(f)
})

# 4. newdata carrying the same levels in a different order
cat("\n-- newdata whose factor has the levels in another order\n")
f <- frm_g("gf")
nd <- d[c(2L, 1L), , drop = FALSE]
nd$gf <- factor(as.character(nd$gf), levels = c("b", "a"))
say("levels reversed on newdata", ndt_time(f, nd))
nd2 <- d[c(2L, 1L), , drop = FALSE]
say("levels as fitted", ndt_time(f, nd2))
cat("   group floors:",
    paste(sprintf("%.6f", tapply(d$rt, d$gf, min)), collapse = " "),
    "\n")
