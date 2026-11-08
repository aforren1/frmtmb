# Review of lane loglik, short check of fix round 1: df = NA under a
# density prior and its readers (M2), uniform priors as bounds (m3), the
# resolved_prior slot (m8), and the restart and allfit labels (m6).
# Lane library first. Seeds 701 to 706.
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), format(packageVersion("frmtmb")),
    "\n")
try_v <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(paste0("ERROR <", class(e)[1], "> ",
                       substr(gsub("\\s+", " ", conditionMessage(e)), 1, 150)),
                class = "err")
    }),
    warning = function(x) {
      w <<- c(w, substr(gsub("\\s+", " ", conditionMessage(x)), 1, 90))
      invokeRestart("muffleWarning")
    },
    message = function(x) invokeRestart("muffleMessage"))
  out <- if (inherits(v, "err")) unclass(v) else
    paste(utils::capture.output(print(v)), collapse = " | ")
  paste0(out, if (length(w)) paste0("  [warnings: ",
                                    paste(w, collapse = " // "), "]"))
}
say <- function(lbl, expr) cat(sprintf("%-40s %s\n", lbl, try_v(expr)))

cat("\n## 1. df = NA and its readers\n")
set.seed(701)
d <- data.frame(x = rnorm(120), z = rnorm(120),
                g = factor(rep(1:20, each = 6)))
d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] + rnorm(120)
pr <- set_prior("normal(0, 0.1)", class = "b")
fp <- frm(y ~ x + (1 | g), data = d, prior = pr)
fp2 <- frm(y ~ x + z + (1 | g), data = d, prior = pr)
f0 <- frm(y ~ x + (1 | g), data = d)
f02 <- frm(y ~ x + z + (1 | g), data = d)
lmf <- lm(y ~ x, data = d)
say("print(logLik(fp))", logLik(fp))
say("attr df_why", attr(logLik(fp), "df_why"))
say("print(logLik(f0))", logLik(f0))
say("AIC(lm, fp)", AIC(lmf, fp))
say("BIC(lm, fp)", BIC(lmf, fp))
say("AIC(logLik(fp))", AIC(logLik(fp)))
say("AIC(fp)", AIC(fp))
say("anova(fp, fp2)", anova(fp, fp2))
say("anova(f0, fp)", anova(f0, fp))
say("anova(fp)", anova(fp))
say("drop1(fp)", drop1(fp))
say("step(fp)", stats::step(fp, trace = 0))
say("lmtest::lrtest(fp, fp2)", lmtest::lrtest(fp, fp2))
say("lmtest::lrtest(f0, fp2)", lmtest::lrtest(f0, fp2))
say("lmtest::lrtest(f0, f02) control", lmtest::lrtest(f0, f02))
say("insight::get_loglikelihood(fp)", insight::get_loglikelihood(fp))
say("insight::get_df(fp, 'model')", insight::get_df(fp, type = "model"))
say("bbmle::AICtab(f0, fp)", bbmle::AICtab(f0, fp))
say("bbmle::AICctab(f0, fp)", bbmle::AICctab(f0, fp))
say("bbmle::AICtab(f0, f02) control", bbmle::AICtab(f0, f02))
say("deviance(fp)", deviance(fp))
say("df.residual(fp)", df.residual(fp))
say("summary(fp) Method",
    grep("Method", utils::capture.output(print(summary(fp))), value = TRUE))
say("print(fp) runs", length(utils::capture.output(print(fp))))
af <- frm_allfit(fp, optimizers = list(nlminb = "nlminb", optim = "optim"))
cat("print(frm_allfit(fp)):\n")
cat(paste0("   ", utils::capture.output(print(af))), sep = "\n")
set.seed(702)
imps <- lapply(1:3, function(i) {
  dd <- data.frame(x = rnorm(60), z = rnorm(60))
  dd$y <- 1 + 0.5 * dd$x + 0.2 * dd$z + rnorm(60)
  dd
})
m1 <- frm_multiple(bf(y ~ x + z), data = imps, prior = pr)
m0 <- frm_multiple(bf(y ~ x), data = imps, prior = pr)
n1 <- frm_multiple(bf(y ~ x + z), data = imps)
n0 <- frm_multiple(bf(y ~ x), data = imps)
say("anova(m1, m0, D1) prior", anova(m1, m0, method = "D1"))
say("anova(n1, n0, D1) no prior", anova(n1, n0, method = "D1"))
say("anova(m1, m0, D2 wald) prior", anova(m1, m0, method = "D2", use = "wald"))
say("anova(m1, m0, D3) prior", anova(m1, m0))
say("anova(m1, constraint = 'z', D1)", anova(m1, constraint = "z",
                                              method = "D1"))
say("D1: pooled z est / se, prior", {
  s <- summary(m1)
  s$fixed[grep("z", rownames(s$fixed)), 1:2, drop = FALSE]
})

cat("\n## 2. uniform priors\n")
set.seed(703)
d3 <- data.frame(x = rnorm(120), g = factor(rep(1:20, each = 6)))
d3$y <- 1 + 0.5 * d3$x + rnorm(20, 0, 0.7)[d3$g] + rnorm(120)
d3$o <- factor(cut(0.8 * d3$x + rlogis(120), c(-Inf, -1, 0, 1, Inf),
                   labels = FALSE), ordered = TRUE)
fit_row <- function(lbl, ..., ref = NULL) {
  f <- tryCatch(withCallingHandlers(frm(...),
                  warning = function(w) invokeRestart("muffleWarning"),
                  message = function(m) invokeRestart("muffleMessage")),
                error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("%-40s FIT ERROR %s\n", lbl, substr(conditionMessage(f), 1, 120)))
    return(invisible(NULL))
  }
  cat(sprintf("%-40s b_x %s  df %s  AIC %s  log_prior %s  method %s\n", lbl,
              if ("x" %in% names(f$estimates$beta))
                format(f$estimates$beta[["x"]], digits = 9) else "-",
              format(attr(logLik(f), "df")),
              substr(try_v(AIC(f)), 1, 50),
              format(attr(logLik_penalized(f), "log_prior"), digits = 6),
              summary(f)$algorithm))
  invisible(f)
}
fit_row("no prior", y ~ x + (1 | g), data = d3)
fit_row("uniform(-5, 5) b x", y ~ x + (1 | g), data = d3,
        prior = set_prior("uniform(-5, 5)", class = "b", coef = "x"))
fit_row("uniform(-5,5) b x, lb/ub", y ~ x + (1 | g), data = d3,
        prior = set_prior("uniform(-5, 5)", class = "b", coef = "x",
                          lb = -5, ub = 5))
fit_row("BINDING uniform(0, 0.2) b x, lb/ub", y ~ x + (1 | g), data = d3,
        prior = set_prior("uniform(0, 0.2)", class = "b", coef = "x",
                          lb = 0, ub = 0.2))
fit_row("BINDING bounds only lb 0 ub 0.2", y ~ x + (1 | g), data = d3,
        prior = set_prior("", class = "b", coef = "x", lb = 0, ub = 0.2))
fit_row("BINDING uniform(0, 0.2), no lb/ub", y ~ x + (1 | g), data = d3,
        prior = set_prior("uniform(0, 0.2)", class = "b", coef = "x"))
cat("   (ML AIC for reference: the no-prior row)\n")
fit_row("uniform(0, 5) sd", y ~ x + (1 | g), data = d3,
        prior = set_prior("uniform(0, 5)", class = "sd"))
fit_row("uniform(0, 5) sigma", y ~ x + (1 | g), data = d3,
        prior = set_prior("uniform(0, 5)", class = "sigma"))
fit_row("uniform(-3, 3) Intercept dpar sigma", bf(y ~ x + (1 | g), sigma ~ 1),
        data = d3,
        prior = set_prior("uniform(-3, 3)", class = "Intercept",
                          dpar = "sigma"))
fit_row("uniform(-3, 3) b dpar sigma", bf(y ~ x, sigma ~ x), data = d3,
        prior = set_prior("uniform(-3, 3)", class = "b", dpar = "sigma"))
fit_row("uniform(-5, 5) Intercept (centered)", y ~ x + (1 | g), data = d3,
        prior = set_prior("uniform(-5, 5)", class = "Intercept"))
fit_row("uniform(-5, 5) thresholds", o ~ x, data = d3, family = cumulative(),
        prior = set_prior("uniform(-5, 5)", class = "Intercept"))
d3$t <- rep(1:6, 20)
fit_row("uniform(-0.9, 0.9) ar", y ~ x + ar(time = t, gr = g, p = 1),
        data = d3, prior = set_prior("uniform(-0.9, 0.9)", class = "ar"))
fit_row("uniform(0, 3) nlpar k", bf(y ~ a + k * x, a + k ~ 1, nl = TRUE),
        data = d3, prior = set_prior("uniform(0, 3)", nlpar = "k", lb = 0,
                                     ub = 3) +
          set_prior("normal(0, 10)", nlpar = "a"))
fit_row("uniform(0, 3) nlpar k alone", bf(y ~ a + k * x, a + k ~ 1, nl = TRUE),
        data = d3, prior = set_prior("uniform(0, 3)", nlpar = "k", lb = 0,
                                     ub = 3))

fa <- suppressWarnings(frm(y ~ x + ar(time = t, gr = g, p = 1), data = d3,
                           prior = set_prior("normal(0, 0.3)", class = "ar")))
cat("ar prior entry: comp", fa$resolved_prior[[1]]$comp, " scale",
    format(fa$resolved_prior[[1]]$scale), " kind",
    fa$resolved_prior[[1]]$dist$kind, "\n")

cat("\n## 5. the resolved_prior slot\n")
lp_of <- function(f, ent) -as.numeric(frmtmb:::neg_log_prior_fn(ent)(f$estimates))
fresh <- function(f) frmtmb:::resolve_prior_input(
  list(frame = f$frame, spec = f$spec), f$prior)$entries
chk <- function(lbl, f) {
  a <- lp_of(f, f$resolved_prior)
  b <- lp_of(f, fresh(f))
  cat(sprintf("%-40s slot %s  slot log prior %.10f  fresh %.10f  same %s\n",
              lbl, "resolved_prior" %in% names(f), a, b, identical(a, b)))
}
# a data-dependent offset: the centered Intercept's prior reads means_X
pc <- set_prior("normal(1, 0.5)", class = "Intercept") +
  set_prior("normal(0, 0.2)", class = "b") +
  set_prior("normal(0, 0.5)", class = "sd")
fc <- frm(y ~ x + (1 | g), data = d, prior = pc)
chk("fit", fc)
d2 <- d
d2$x <- d2$x + 3
chk("update(data = shifted x)", update(fc, data = d2))
chk("update(. ~ . + z)", update(fc, . ~ . + z))
chk("update(prior = new)", update(fc, prior = set_prior("normal(0, 1)",
                                                        class = "b")))
rf <- try(lme4::refit(fc, newresp = d$y + rnorm(nrow(d), 0, 0.1)))
if (!inherits(rf, "try-error")) chk("lme4::refit(newresp)", rf)
mc <- frm_multiple(bf(y ~ x + z), data = lapply(imps, function(dd) {
  dd$x <- dd$x + runif(1, 0, 3); dd }), prior = pc[1:2])
for (j in seq_along(mc$fits)) chk(paste("frm_multiple fit", j), mc$fits[[j]])
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
            paste(grep("^(prior|res)", names(fc), value = TRUE), collapse = ",")))
cat(sprintf("no-prior fit: slot present %s, value NULL %s, density %s\n",
            "resolved_prior" %in% names(f0), is.null(f0$resolved_prior),
            frmtmb:::fit_prior_density(f0)))
# an object without the slot falls back to resolving again
fs <- fc
fs$resolved_prior <- NULL
cat(sprintf("slot removed: logLik identical %s\n",
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
    message = function(m) { msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage") },
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning") })
  cat(sprintf("seed %d logLik %.6f  logLik_penalized %.6f\n", seed,
              as.numeric(logLik(f)), as.numeric(logLik_penalized(f))))
  for (m in grep("log-likelihood", msgs, value = TRUE)) {
    cat("   ", gsub("\\s+", " ", m), "\n")
  }
}
cat("done\n")
