# Review of lane loglik, short check, part 2: sections 5 and 6 of
# review-loglik-short.R (which stopped on its own frm_multiple() prior
# bug), the AR chart's log-Jacobian value, and insight on a no-prior fit.
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), format(packageVersion("frmtmb")),
    "\n")
set.seed(701)
d <- data.frame(x = rnorm(120), z = rnorm(120),
                g = factor(rep(1:20, each = 6)))
d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] + rnorm(120)
f0 <- frm(y ~ x + (1 | g), data = d)
fp <- frm(y ~ x + (1 | g), data = d,
          prior = set_prior("normal(0, 0.1)", class = "b"))
set.seed(702)
imps <- lapply(1:3, function(i) {
  dd <- data.frame(x = rnorm(60), z = rnorm(60))
  dd$y <- 1 + 0.5 * dd$x + 0.2 * dd$z + rnorm(60)
  dd
})

cat("\n## 5. the resolved_prior slot (continued)\n")
lp_of <- function(f, ent) {
  -as.numeric(frmtmb:::neg_log_prior_fn(ent)(f$estimates))
}
fresh <- function(f) frmtmb:::resolve_prior_input(
  list(frame = f$frame, spec = f$spec), f$prior)$entries
chk <- function(lbl, f) {
  a <- lp_of(f, f$resolved_prior)
  b <- lp_of(f, fresh(f))
  cat(sprintf("%-32s slot log prior %.10f  fresh %.10f  same %s\n",
              lbl, a, b, identical(a, b)))
}
pc2 <- set_prior("normal(1, 0.5)", class = "Intercept") +
  set_prior("normal(0, 0.2)", class = "b")
mc <- frm_multiple(bf(y ~ x + z), data = lapply(seq_along(imps), function(j) {
  dd <- imps[[j]]
  dd$x <- dd$x + j
  dd
}), prior = pc2)
for (j in seq_along(mc$fits)) chk(paste("frm_multiple fit", j), mc$fits[[j]])
fc <- frm(y ~ x + (1 | g), data = d,
          prior = pc2 + set_prior("normal(0, 0.5)", class = "sd"))
tf <- tempfile(fileext = ".rds")
saveRDS(fc, tf)
fr <- readRDS(tf)
cat(sprintf("readRDS: logLik identical %s, logLik_penalized identical %s, df NA %s\n",
            identical(logLik(fr), logLik(fc)),
            identical(logLik_penalized(fr), logLik_penalized(fc)),
            is.na(attr(logLik(fr), "df"))))
unlink(tf)
cat(sprintf("partial match: identical(fc$prior, fc[['prior']]) %s; names starting 'prior' or 'res': %s\n",
            identical(fc$prior, fc[["prior"]]),
            paste(grep("^(prior|res)", names(fc), value = TRUE),
                  collapse = ",")))
cat(sprintf("no-prior fit: slot present %s, value NULL %s, density %s\n",
            "resolved_prior" %in% names(f0), is.null(f0$resolved_prior),
            frmtmb:::fit_prior_density(f0)))
fs <- fc
fs$resolved_prior <- NULL
cat(sprintf("slot removed: names has slot %s, logLik identical %s\n",
            "resolved_prior" %in% names(fs),
            identical(as.numeric(logLik(fs)), as.numeric(logLik(fc)))))

cat("\n## 6. labels\n")
for (seed in 501:502) {
  set.seed(seed)
  dn <- data.frame(x = runif(80, 0, 3))
  dn$y <- 0.8 * dn$x + rnorm(80, 0, 0.3)
  msgs <- character(0)
  f <- withCallingHandlers(
    frm(bf(y ~ a * b * x, a + b ~ 1, nl = TRUE), data = dn,
        prior = set_prior("normal(0, 1)", nlpar = "a") +
          set_prior("normal(0, 1)", nlpar = "b")),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  cat(sprintf("seed %d logLik %.6f  logLik_penalized %.6f\n", seed,
              as.numeric(logLik(f)), as.numeric(logLik_penalized(f))))
  for (m in grep("log-likelihood", msgs, value = TRUE)) {
    cat("   ", gsub("\\s+", " ", m), "\n")
  }
}

cat("\n## extra\n")
set.seed(703)
d3 <- data.frame(x = rnorm(120), g = factor(rep(1:20, each = 6)),
                 t = rep(1:6, 20))
d3$y <- 1 + 0.5 * d3$x + rnorm(20, 0, 0.7)[d3$g] + rnorm(120)
fa <- suppressWarnings(frm(y ~ x + ar(time = t, gr = g, p = 1), data = d3,
                           prior = set_prior("normal(0, 0.3)", class = "ar")))
e <- fa$resolved_prior[[1]]
cat("ar thetaac", fa$estimates$thetaac, " trans", format(e$dist$trans),
    " logjac", frmtmb:::ac_trans_logjac(fa$estimates$thetaac[e$idx],
                                         e$dist$trans), "\n")
print(insight::get_loglikelihood(f0))
print(insight::get_loglikelihood(fp))
cat("insight::get_df(f0, model):", insight::get_df(f0, type = "model"),
    " attr(logLik(f0), df):", attr(logLik(f0), "df"), "\n")
cat("done\n")
