## REVIEW claim 5: defect B. Look for ANY route besides the ordinal
## threshold vector whose entry covers more than one parameter, and check
## what each one did before and does now: class "cor" (an LKJ density on a
## block's correlation thetas), class "ar" at order > 1 (thetaac), class
## "rescor" (thetar). Also confirm the FIT path and the prior reports are
## untouched: prior_summary(), default_prior() and the prior objective.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, " ver =",
    as.character(utils::packageVersion("frmtmb")), "\n\n")
options(width = 140)

say <- function(tag, expr) {
  cat("==", tag, "\n")
  out <- tryCatch(expr, error = function(e) e, warning = function(w) w)
  if (inherits(out, "condition")) {
    cat("   ", class(out)[1L], ":",
        gsub("\n", " ", conditionMessage(out)), "\n")
  } else {
    cat("    OK\n")
    print(out)
  }
  cat("\n")
  invisible(out)
}

## ---- 1. class "cor": an LKJ density over three correlation thetas ---
set.seed(801)
n <- 90
d <- data.frame(g = factor(rep(1:15, each = 6)), x = rnorm(n), z = rnorm(n))
d$y <- rnorm(n)
say("frm_simulate, (1 + x + z | g), prior on class cor (idx length 3)",
    {
      s <- frm_simulate(bf(y ~ x + z + (1 + x + z | g)), family = gaussian(),
                        data = d,
                        newparams = list(Intercept = 0, x = 1, z = 1,
                                         sigma = 1,
                                         `sd_g__Intercept` = 1,
                                         `sd_g__x` = 1, `sd_g__z` = 1),
                        prior = set_prior("lkj(2)", class = "cor"),
                        nsim = 2, seed = 810)
      attr(s, "pars")
    })

## ---- 2. class "ar" at order 2: two thetaac parameters ---------------
set.seed(802)
d2 <- data.frame(t = 1:60, gg = factor(rep(1, 60)))
d2$x <- rnorm(60)
d2$y <- rnorm(60)
say("frm_simulate, ar(p = 2), prior on class ar (idx length 2)",
    {
      s <- frm_simulate(bf(y ~ x + ar(time = t, gr = gg, p = 2)),
                        family = gaussian(), data = d2,
                        newparams = list(Intercept = 0, x = 1, sigma = 1),
                        prior = set_prior("normal(0, 0.3)", class = "ar"),
                        nsim = 2, seed = 811)
      attr(s, "pars")
    })

## ---- 3. class "rescor": is a multivariate frm_simulate refused? -----
say("frm_simulate on a multivariate model (rescor route)",
    {
      dm <- data.frame(x = rnorm(40), y1 = rnorm(40), y2 = rnorm(40))
      frm_simulate(bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE),
                   family = gaussian(), data = dm,
                   newparams = list(sigma = 1),
                   prior = set_prior("lkj(2)", class = "rescor"),
                   nsim = 2, seed = 812)
    })

## ---- 4. the ordinal threshold vector, the lane's own case -----------
set.seed(803)
dd <- data.frame(x = rnorm(40), y = rep(1:4, 10))
for (fn in c("cumulative", "sratio", "cratio", "acat")) {
  fam <- switch(fn, cumulative = cumulative(), sratio = sratio(),
                cratio = cratio(), acat = acat())
  say(paste0("frm_simulate, ", fn, ", prior on class Intercept (K = 4)"),
      {
        s <- frm_simulate(bf(y ~ x), family = fam, data = dd,
                          prior = set_prior("normal(0, 2)",
                                            class = "Intercept") +
                            set_prior("normal(0, 1)", class = "b"),
                          nsim = 2, seed = 702)
        attr(s, "pars")
      })
}
## two categories is ONE threshold: the draw must still happen
set.seed(804)
dd2 <- data.frame(x = rnorm(40), y = rep(1:2, 20))
for (fn in c("cratio", "acat", "sratio")) {
  fam <- switch(fn, cratio = cratio(), acat = acat(), sratio = sratio())
  say(paste0("frm_simulate, ", fn, ", K = 2 (one threshold)"),
      {
        s <- frm_simulate(bf(y ~ x), family = fam, data = dd2,
                          prior = set_prior("normal(0, 2)",
                                            class = "Intercept") +
                            set_prior("normal(0, 1)", class = "b"),
                          nsim = 2, seed = 704)
        attr(s, "pars")
      })
}

## ---- 5. the FIT path must be untouched -----------------------------
cat("#### the prior objective and the prior reports\n")
for (fn in c("cumulative", "sratio")) {
  fam <- switch(fn, cumulative = cumulative(), sratio = sratio())
  fit <- suppressWarnings(
    frm(bf(y ~ x), family = fam, data = dd,
        prior = set_prior("normal(0, 2)", class = "Intercept") +
          set_prior("normal(0, 1)", class = "b")))
  cat("==", fn, ": logLik =", sprintf("%.14f", logLik(fit)),
      " objective =", sprintf("%.14f", fit$opt$objective), "\n")
  cat("   coefs:",
      paste(sprintf("%s=%.14g", names(frmtmb:::get_coef.frmtmb_fit(fit)),
                    frmtmb:::get_coef.frmtmb_fit(fit)),
            collapse = "  "), "\n")
  ps <- prior_summary(fit)
  cat("   prior_summary rows:", nrow(as.data.frame(ps)), "\n")
  print(as.data.frame(ps))
}
dp <- default_prior(bf(y ~ x), family = cumulative(), data = dd)
cat("== default_prior rows:", nrow(as.data.frame(dp)), "\n")
print(as.data.frame(dp))
cat("DONE rev-06\n")
