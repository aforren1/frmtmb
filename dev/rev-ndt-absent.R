# rev-ndt round 2, attacks 1b and 2: the absent group, and the
# constructed collision used against the fit.
#
# `match()`/`[` on a name that is not in the table returns NA, and an NA
# bound flowing into a density is silent, so every path that can reach
# the scaler with a group the fit never saw is exercised here. The
# collision found by dev/rev-ndt-collide2.R is then used twice: once
# with both labels in ONE column, where ddm_coerce_ndt_group() can see
# them, and once with the second label only on newdata, where the code
# says it cannot.
#
# Seed 4242. Worktree build at rev-ndt-lib2.

.libPaths(c("C:/Users/adf44/source/r/rev-ndt-lib2",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

say <- function(lbl, want, expr) {
  r <- tryCatch({
    v <- force(expr)
    list(ok = TRUE, txt = paste(sprintf("%.6f",
      utils::head(as.numeric(v), 3)), collapse = " "))
  }, error = function(e)
    list(ok = FALSE, txt = substr(conditionMessage(e), 1, 88)))
  got <- if (r$ok) "ran" else "refused"
  flag <- if (identical(got, want)) "   " else "** "
  cat(sprintf("%s%-44s want=%-7s got=%-7s %s\n", flag, lbl, want, got,
              r$txt))
  invisible(r)
}

set.seed(4242)
d <- ddm_simulate(600, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
d$g <- factor(rep(c("s1", "s2", "s3"), length.out = nrow(d)))
form <- bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
           bias = 0.5)
fit <- frm(form, family = wiener(), data = d)
fl <- tapply(d$rt, d$g, min)
cat("group floors:",
    paste(names(fl), sprintf("%.6f", fl), sep = "=", collapse = "  "),
    "\n\n")

one <- d[match(levels(d$g), as.character(d$g)), , drop = FALSE]
unseen <- one[1L, , drop = FALSE]
unseen$g <- factor("s9")
gone <- one
gone$g <- NULL

cat("-- 2. the absent group, on every path that can reach the scaler\n")
say("ndt_time(), unseen group", "refused", ndt_time(fit, unseen))
say("ndt_time(), column gone", "refused", ndt_time(fit, gone))
say("predict(dpar=ndt), unseen group", "refused",
    predict(fit, newdata = unseen, dpar = "ndt", type = "response"))
say("predict(type=response), unseen group", "refused",
    predict(fit, newdata = unseen, type = "response"))
say("predict(type=response), column gone", "refused",
    predict(fit, newdata = gone, type = "response"))
say("fitted mean on a good newdata row", "ran",
    predict(fit, newdata = one, type = "response"))
say("frm() refit on rows naming an unseen group", "refused", {
  d2 <- rbind(d, transform(d[1L, , drop = FALSE], g = "s9"))
  d2$g <- factor(as.character(d2$g))
  frm(form, family = wiener(), data = d2,
      start = NULL)
})
say("residuals(newdata=) with an unseen group", "refused",
    residuals(fit, newdata = unseen, type = "response"))
say("simulate() on the fit", "ran",
    { set.seed(3); unlist(stats::simulate(fit, nsim = 1L)[[1L]][1:2]) })
say("frm_simulate() with an unseen group", "refused", {
  d3 <- one
  d3$g <- factor(c("s9", "s2", "s3"))
  frm_simulate(form, data = d3, family = family(fit),
               newparams = list(beta = 1, betad = c(log(1.4), 0, 0)),
               nsim = 1)
})
say("influence() on the fit", "ran",
    utils::head(as.numeric(unlist(frmtmb::influence(fit,
      group = "g")[1L, , drop = TRUE])), 3))

# --------------------------------------------------------------------
cat("\n-- 1b. the constructed collision, both labels in ONE column\n")
cl <- readRDS(paste0("C:/Users/adf44/source/r/frmtmb-wt-ndt/dev/",
                     "rev-ndt-collide2.rds"))
pkg_code <- get("ddm_label_code", envir = asNamespace("frmtmb.eam"))
cat("   label A code points: ", paste(cl$ca, collapse = " "), "\n",
    sep = "")
cat("   label B code points: ", paste(cl$cb, collapse = " "), "\n",
    sep = "")
cat("   codes: ", sprintf("%.0f", pkg_code(cl$a)), " and ",
    sprintf("%.0f", pkg_code(cl$b)), " equal=",
    identical(pkg_code(cl$a), pkg_code(cl$b)), "\n", sep = "")

set.seed(4242)
dc <- ddm_simulate(600, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
dc$g <- factor(rep(c(cl$a, cl$b, "s3"), length.out = nrow(dc)))
say("fit with BOTH colliding labels present", "refused",
    frm(form, family = wiener(), data = dc))

cat("\n-- 1c. label A fitted, label B seen only on newdata\n")
set.seed(4242)
da <- ddm_simulate(600, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
da$g <- factor(rep(c(cl$a, "s3"), length.out = nrow(da)))
fa <- frm(form, family = wiener(), data = da)
fla <- tapply(da$rt, da$g, min)
cat("   fitted floors: ", paste(sprintf("%.6f", fla), collapse = "  "),
    "\n", sep = "")
ndA <- da[match(cl$a, as.character(da$g)), , drop = FALSE]
ndB <- ndA
ndB$g <- factor(cl$b)
ndU <- ndA
ndU$g <- factor("s9")
tA <- say("ndt_time() on the FITTED label A", "ran", ndt_time(fa, ndA))
tB <- say("ndt_time() on the COLLIDING label B", "refused",
          ndt_time(fa, ndB))
tU <- say("ndt_time() on an ordinary unseen label", "refused",
          ndt_time(fa, ndU))
if (tB$ok && tA$ok) {
  cat("   label B was scored against label A's bound: ",
      identical(tA$txt, tB$txt), "\n", sep = "")
}
say("predict(type=response) on label B", "refused",
    predict(fa, newdata = ndB, type = "response"))
