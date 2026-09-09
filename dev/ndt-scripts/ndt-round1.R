# ndt lane, round 1: the five blockers and the nits, measured.
#
#   Rscript --vanilla dev/ndt-scripts/ndt-round1.R ref|new
#
# The `ref` arm is the 0.6.0 build at reflib-r2 and skips everything
# that needs ndt_group() or ndt_time(). Seeds are named per block and
# follow the review's: 4242 for the prior and the pinned constant, 808
# for the level ordering.
a <- commandArgs(trailingOnly = TRUE)
arm <- if (length(a)) a[1L] else "new"
.libPaths(if (identical(arm, "ref")) {
  c("C:/Users/adf44/source/r/reflib-r2",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
} else {
  c("C:/Users/adf44/source/r/ndt-lib",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
})
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("arm:", arm, "\n\n")
fmt <- function(x, d = 6) formatC(x, digits = d, format = "g")

## ------------------------------------------------- the scalar path (B2, B3)
set.seed(4242)
d <- ddm_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
d$cond <- factor(rep(0:1, each = 200), labels = c("a", "b"))
cat("min(rt):", fmt(min(d$rt), 8), "\n")

f0 <- bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5)
fit <- frm(f0, family = wiener(), data = d)
cat("ndt link name        :", family(fit)$links$ndt$name, "\n")
cat("predict(ndt, response):",
    fmt(predict(fit, dpar = "ndt", type = "response")[1], 8), "\n")
if (!identical(arm, "ref")) {
  cat("ndt_time()            :", fmt(ndt_time(fit)[1], 8), "\n")
}
cat("logLik                :", fmt(as.numeric(logLik(fit)), 12), "\n")

# B2: a prior written the way a brms user writes it. class = "<dpar>"
# applies only where the dpar has no formula of its own, which is brms's
# rule and frmtmb's, so it gets a model that gives `ndt` no predictor.
fnp <- bf(rt | dec(upper) ~ cond, bs ~ 1, bias = 0.5)
for (loc in c(0.30, 0.10)) {
  fp <- frm(fnp, family = wiener(), data = d,
            prior = prior_string(paste0("normal(", loc, ", 0.01)"),
                                 class = "ndt"))
  cat("prior normal(", loc, ", 0.01) class=ndt -> ndt =",
      fmt(predict(fp, dpar = "ndt", type = "response")[1], 8), "\n")
}
fi <- frm(f0, family = wiener(), data = d,
          prior = prior_string("normal(0, 0.1)",
                               class = "Intercept", dpar = "ndt"))
cat("prior on the Intercept, dpar=ndt -> eta =",
    fmt(unlist(fixef(fi))[["ndt.(Intercept)"]], 8), "\n")

# B3: a pinned constant, with and without max_ndt
fc <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt = 0.2, bias = 0.5),
          family = wiener(max_ndt = 0.25), data = d)
cat("pinned ndt = 0.2 with max_ndt = 0.25: logLik",
    fmt(as.numeric(logLik(fc)), 10), " fitted[1]",
    fmt(fitted(fc)[1], 6), "\n")
cat("pinned ndt = 0.2 with a bare wiener(): ",
    tryCatch({
      frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt = 0.2, bias = 0.5),
          family = wiener(), data = d)
      "NO ERROR"
    }, error = function(e) paste("refused:",
                                 substr(conditionMessage(e), 1, 48))),
    "\n")

if (identical(arm, "ref")) quit(save = "no")

## ---------------------------------------------------------- B1: level order
set.seed(808)
dg <- ddm_simulate(900, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
dg$g <- factor(rep(c("s1", "s2", "s3"), length.out = nrow(dg)))
cat("\ngroup floors:",
    paste(names(tapply(dg$rt, dg$g, min)),
          fmt(as.numeric(tapply(dg$rt, dg$g, min))), sep = "=",
          collapse = "  "), "\n")
fg <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
             bias = 0.5), family = wiener(), data = dg)
bd <- frmtmb::single_response(fg)[["family"]][["ndt_bound"]]
cat("floors captured, keyed by LABEL CODE:", fmt(bd$floors), "\n")
cat("group sizes recorded               :", bd$sizes, "\n")

sub <- dg[dg$g != "s1", ]
grids <- list(
  "subset, levels as fitted" = sub,
  "subset then droplevels()" = droplevels(sub),
  "subset then relevel(s3)"  = transform(sub, g = relevel(g, "s3")),
  "hand-built, reversed"     = data.frame(
    rt = sub$rt, upper = sub$upper,
    g = factor(as.character(sub$g), levels = c("s3", "s2", "s1"))))
ref <- NULL
for (nm in names(grids)) {
  nd <- grids[[nm]]
  one <- nd[match(c("s2", "s3"), as.character(nd$g)), , drop = FALSE]
  v <- as.numeric(ndt_time(fg, newdata = one))
  m <- as.numeric(predict(fg, newdata = one, type = "response"))
  if (is.null(ref)) ref <- v
  cat(sprintf("  %-26s ndt_time %s  mean %s  same as first: %s\n",
              nm, paste(fmt(v), collapse = " "),
              paste(fmt(m), collapse = " "),
              identical(all.equal(v, ref), TRUE)))
}
# a character and an integer column name the same groups, so they must
# give the same fit
for (sp in c("chr", "int", "lgl")) {
  d2 <- dg
  d2$gg <- switch(sp, chr = as.character(d2$g),
                  int = as.integer(d2$g), lgl = d2$g == "s2")
  f2 <- frm(bf(rt | dec(upper) + ndt_group(gg) ~ 1, bs ~ 1, ndt ~ 1,
               bias = 0.5), family = wiener(), data = d2)
  cat("  spelling", sp, "logLik", fmt(as.numeric(logLik(f2)), 10),
      " floors", fmt(sort(frmtmb::single_response(f2)[["family"]][[
        "ndt_bound"]]$floors)), "\n")
}

## ------------------------- B1b: a label the floor table does not hold
## match() returns NA on a miss and an NA bound flowing into the density
## would be the same defect one layer down, so the case is constructed.
cat("\nB1b, a group label the fit never saw:\n")
nd_new <- dg[1:3, ]
nd_new$g <- factor(rep("s9", 3), levels = c(levels(dg$g), "s9"))
for (what in c("ndt_time", "predict(response)", "frm refit")) {
  msg <- tryCatch({
    switch(what,
      "ndt_time" = ndt_time(fg, newdata = nd_new),
      "predict(response)" = predict(fg, newdata = nd_new,
                                    type = "response"),
      "frm refit" = frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1,
                           ndt ~ 1, bias = 0.5),
                        family = frmtmb::single_response(fg)[["family"]],
                        data = nd_new))
    "NO ERROR (fails open)"
  }, error = function(e) paste("refused:",
                               substr(conditionMessage(e), 1, 58)))
  cat(sprintf("  %-18s %s\n", what, msg))
}

## -------------------------------------------------- N0: the bound is KEPT
cat("\nN0, refinalize without the fastest row:\n")
set.seed(4242)
dn <- ddm_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
dn$w <- rdm_simulate(400, v = c(2.5, 1.5), A = 0.5, k = 0.5,
                     ndt = 0.2)$choice
drop <- -which.min(dn$rt)
cases <- list(
  wiener = list(f = wiener(),
                at = list(dec = as.numeric(dn$upper))),
  lba = list(f = lba(2), at = list(vint1 = as.numeric(dn$w))),
  rdm = list(f = rdm(2), at = list(vint1 = as.numeric(dn$w))),
  wiener_gng = list(f = wiener_gng(deadline = 5),
                    at = list(dec = rep(1, nrow(dn)))))
for (nm in names(cases)) {
  fam <- cases[[nm]]$f
  at <- cases[[nm]]$at
  fin <- fam[["family_finalize"]](fam, dn$rt, at)
  at2 <- lapply(at, function(v) v[drop])
  fin2 <- fin[["family_finalize"]](fin, dn$rt[drop], at2)
  cat(sprintf("  %-11s all rows %s  refinalized %s  KEPT: %s\n", nm,
              fmt(fin[["ndt_bound"]][["ub"]], 8),
              fmt(fin2[["ndt_bound"]][["ub"]], 8),
              identical(fin[["ndt_bound"]][["ub"]],
                        fin2[["ndt_bound"]][["ub"]])))
}

## ------------------------------------------------------- N2 and N4 and N5
cat("\nN2, ndt_time() on newdata with no group column: ",
    tryCatch({
      ndt_time(fg, newdata = data.frame(rt = 0.5, upper = 1))
      "NO ERROR"
    }, error = function(e) paste("refused:",
                                 substr(conditionMessage(e), 1, 60))),
    "\n")
cat("N4, a group of one trial, its recorded size: ")
d1 <- dg
d1$g <- factor(c("solo", as.character(dg$g[-1])))
f1 <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
             bias = 0.5), family = wiener(), data = d1)
b1 <- frmtmb::single_response(f1)[["family"]][["ndt_bound"]]
cat(min(b1$sizes), "\n")
cat("N5, an NA group through ndt_time(newdata =): ",
    tryCatch({
      nd <- dg[1:2, ]
      nd$g[1] <- NA
      ndt_time(fg, newdata = nd)
      "NO ERROR"
    }, error = function(e) paste("refused:",
                                 substr(conditionMessage(e), 1, 45))),
    "\n")
