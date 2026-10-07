# Reviewer of lane setier, re-check of item 8 (the boundary message in
# place of the convergence warning at code != 0). A fit that stopped
# short and is NOT at a boundary must still warn.
#   Rscript dev/setier-rev2-item8.R <lib>
# Cases (seeds 1..10 each):
#   id_iter   y ~ x + (1 | g), group sd 0.8 (identified), eval.max 8
#   id_glm    bernoulli y ~ x + z, eval.max 6
#   bnd_iter  y ~ x + (1 | g), no group variance, eval.max 8
#   bnd_rs    dev/setier-singular.R's rs20, default budget
#   nl_iter   y ~ a * exp(b * x), eval.max 6 (identified)
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), "\n")
short <- function(k) frmtmb_control(optCtrl = list(eval.max = k,
                                                   iter.max = k),
                                    restarts = 0)
for (cs in c("id_iter", "id_glm", "bnd_iter", "bnd_rs", "nl_iter")) {
  tab <- character()
  for (s in 1:10) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120),
                    z = rnorm(120))
    w <- character(); m <- character()
    f <- withCallingHandlers(switch(cs,
      id_iter = {
        d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.8)[d$g] + rnorm(120)
        frm(y ~ x + (1 | g), data = d, control = short(8))
      },
      id_glm = {
        d$y <- rbinom(120, 1, plogis(0.3 + d$x))
        frm(y ~ x + z, family = bernoulli(), data = d, control = short(6))
      },
      bnd_iter = {
        d$y <- 1 + 0.5 * d$x + rnorm(120)
        frm(y ~ x + (1 | g), data = d, control = short(8))
      },
      bnd_rs = {
        d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] + rnorm(120)
        frm(y ~ x + (1 + x | g), data = d)
      },
      nl_iter = {
        d$y <- 2 * exp(0.5 * d$x) + rnorm(120, 0, 0.3)
        frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = d,
            control = short(6), start = list(beta = c(1, 0.1)))
      }),
      warning = function(x) {w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")},
      message = function(x) {m <<- c(m, conditionMessage(x))
        invokeRestart("muffleMessage")})
    g <- max(abs(f$obj$gr(f$opt$par) * (f$par_units %||% 1)))
    lost <- ns$sdr_of(f)$se_lost
    conv <- any(grepl("^Optimizer did not report|gradient", w))
    bm <- any(grepl("^Boundary", m))
    tab <- c(tab, sprintf("code%d %s%s", f$opt$convergence,
                          if (conv) "CONVWARN" else "noconv",
                          if (bm) "+BOUNDARY" else ""))
    if (f$opt$convergence != 0 && !conv) {
      cat(sprintf("   %s seed %d: code %d silent on convergence; max grad %.2g; lost %s; boundary msg %s\n",
                  cs, s, f$opt$convergence, g,
                  paste(names(lost), lost, sep = ":", collapse = ","), bm))
    }
  }
  tb <- table(tab)
  cat(sprintf("%-9s %s\n", cs, paste(names(tb), tb, collapse = " | ")))
}
