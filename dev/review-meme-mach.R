# Review of lane meme: the 0.77.0 machinery on me() models.
#   Rscript dev/review-meme-mach.R > dev/review-meme-out/mach.txt
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-meme/dev/review-meme-shapes.R")
d <- rv_data()
try_ <- function(lab, expr) {
  r <- tryCatch({
    v <- withCallingHandlers(expr, message = function(m) {
      cat("  [message]", conditionMessage(m))
      invokeRestart("muffleMessage")
    }, warning = function(w) {
      cat("  [warning]", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
    if (is.null(v)) "ok (NULL)" else v
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("--", lab, "\n")
  if (is.character(r) && length(r) == 1) cat("  ", r, "\n") else print(r)
  invisible(r)
}
ent <- function(fit, pr) {
  ri <- frmtmb:::resolve_priorlist(fit, pr)
  vapply(ri$entries, function(e) {
    paste0(e$comp, "[", paste(e$idx, collapse = ","), "]:", e$dist$kind)
  }, "")
}
f_me <- bf(y ~ me(x, sx) + me(z, sz))
f_mere <- bf(y ~ me(x, sx) + me(z, sz) + w + (1 + w | g))
des_me <- frmtmb:::prior_design(f_me, d, NULL, list())
des_mere <- frmtmb:::prior_design(f_mere, d, NULL, list())

cat("\n#### 1. one slot per parameter\n")
try_("cor lkj(2) on a me()-only model",
     ent(des_me, set_prior("lkj(2)", class = "cor")))
try_("theta normal(0,1) on a me()-only model",
     ent(des_me, set_prior("normal(0, 1)", class = "theta")))
try_("sd on a me()-only model",
     ent(des_me, set_prior("exponential(1)", class = "sd")))
try_("cov_class_coefs rows of a me() + RE model",
     frmtmb:::cov_class_coefs(des_mere$spec, des_mere$frame)[, c("comp", "idx",
                                                              "class")])
try_("cor lkj(2) + corme lkj(3) + sd + sdme on me() + RE",
     ent(des_mere, set_prior("lkj(2)", class = "cor") +
           set_prior("lkj(3)", class = "corme") +
           set_prior("exponential(1)", class = "sd") +
           set_prior("exponential(2)", class = "sdme")))
try_("theta coef = thetame_1",
     ent(des_me, set_prior("normal(0, 1)", class = "theta",
                           coef = "thetame_1")))
try_("two rows on one slot (sdme class twice)",
     ent(des_me, set_prior("exponential(1)", class = "sdme") +
           set_prior("exponential(2)", class = "sdme")))
try_("sdme coef row twice",
     ent(des_me, set_prior("exponential(1)", class = "sdme", coef = "mex") +
           set_prior("exponential(2)", class = "sdme", coef = "mex")))
try_("corme class twice",
     ent(des_me, set_prior("lkj(2)", class = "corme") +
           set_prior("lkj(3)", class = "corme")))

cat("\n#### 2. route mark\n")
tf <- default_prior(f_me, data = d)
ts <- default_prior(f_me, data = d, route = "sample")
cat("attr route fit:", attr(tf, "route"), " sample:", attr(ts, "route"),
    " subset:", attr(tf[tf$class %in% c("sdme", "corme"), ], "route"),
    " c():", attr(c(tf[tf$class == "sdme", ],
                    set_prior("lkj(2)", class = "corme")), "route"), "\n")
try_("frm_sample(prior = fit-route table)",
     {frm_sample(f_me, data = d, family = gaussian(), prior = tf,
                 chains = 1, iter = 50, refresh = 0); "sampled"})
try_("frm_sample(prior = subset of fit-route me rows)",
     {frm_sample(f_me, data = d, family = gaussian(),
                 prior = tf[tf$class %in% c("sdme", "corme"), ],
                 chains = 1, iter = 50, refresh = 0); "sampled"})
try_("frm(prior = sample-route table): logLik vs hand-written",
     {a <- frm(f_me, data = d, prior = ts)
      b <- frm(f_me, data = d, prior = set_prior("lkj(1)", class = "corme") +
                 set_prior("student_t(3, 0.9, 2.5)", class = "Intercept") +
                 set_prior("student_t(3, 0, 2.5)", class = "sigma"))
      c(table = as.numeric(logLik(a)), hand = as.numeric(logLik(b)))})
try_("frm_sample(prior = sample-route table) entries",
     {rs <- frmtmb.sample:::sample_resolve_priors(
       frm(f_me, data = d, dry_run = "objective"), ts)
      vapply(rs$ri$entries, function(e) paste0(e$comp, ":", e$dist$kind), "")})

cat("\n#### 3. sample_prior = 'only' with a coef-only meanme prior\n")
pb <- set_prior("normal(0, 1)", class = "b") +
  set_prior("exponential(1)", class = "sdme") +
  set_prior("lkj(2)", class = "corme")
try_("only: meanme coef = mez alone",
     {frm_sample(f_me, data = d, family = gaussian(), prior = pb +
                   set_prior("normal(0, 1)", class = "meanme", coef = "mez"),
                 sample_prior = "only", chains = 1, iter = 50, refresh = 0)
      "sampled"})
try_("only: corme flat (prior = '' on corme)",
     {frm_sample(f_me, data = d, family = gaussian(), prior =
                   set_prior("normal(0, 1)", class = "b") +
                   set_prior("normal(0, 1)", class = "meanme") +
                   set_prior("exponential(1)", class = "sdme") +
                   set_prior("", class = "corme"),
                 sample_prior = "only", chains = 1, iter = 50, refresh = 0)
      "sampled"})

cat("\n#### 4. prior draws with a coef bound (frmtmb accepts, brms refuses)\n")
ds <- try_("yes: sdme coef = mex exponential(1), ub = 0.5",
           frm_sample(f_me, data = d, family = gaussian(),
                      prior = set_prior("exponential(1)", class = "sdme",
                                        coef = "mex", ub = 0.5) +
                        set_prior("exponential(1)", class = "sdme"),
                      sample_prior = "yes", chains = 2, iter = 1000,
                      seed = 3, refresh = 0))
if (inherits(ds, "frmtmb_draws")) {
  pd <- prior_draws(ds)
  cat("  prior draw names:", names(pd), "\n")
  m <- as.matrix(ds)
  cat(sprintf("  posterior max sdme_mex %.3f (bound 0.5)\n",
              max(m[, "sdme_mex"])))
  for (nm in grep("sdme", names(pd), value = TRUE)) {
    cat(sprintf("  prior %s: max %.3f  mean %.3f\n", nm, max(pd[[nm]]),
                mean(pd[[nm]])))
  }
  print(as.data.frame(prior_summary(ds)))
  cat("\n#### 5. hypothesis() Evid.Ratio on meanme\n")
}
ds2 <- frm_sample(f_me, data = d, family = gaussian(),
                  prior = set_prior("normal(0, 1)", class = "meanme") +
                    set_prior("exponential(1)", class = "sdme"),
                  sample_prior = "yes", chains = 2, iter = 1000, seed = 3,
                  refresh = 0)
cat("  prior draw names:", names(prior_draws(ds2)), "\n")
try_("prior_draws(variable = 'meanme_mex')",
     names(prior_draws(ds2, variable = "meanme_mex")))
try_("hypothesis meanme_mex = 0",
     hypothesis(ds2, "meanme_mex = 0")$hypothesis)
try_("hypothesis meanme_mex > 0",
     hypothesis(ds2, "meanme_mex > 0")$hypothesis)
try_("hypothesis sdme_mex = 1",
     hypothesis(ds2, "sdme_mex = 1")$hypothesis)

cat("\n#### 6. update()\n")
try_("update(draws, prior = corme lkj(4))",
     {u <- update(ds2, prior = set_prior("lkj(4)", class = "corme"),
                  refresh = 0)
      as.data.frame(prior_summary(u))[, c("prior", "class", "coef")]})
try_("update(draws, formula = . ~ . - me(z, sz))",
     {u <- update(ds2, . ~ . - me(z, sz), refresh = 0)
      list(prior = as.data.frame(prior_summary(u))[, c("prior", "class")],
           pd = names(prior_draws(u)))})
try_("update(frm fit, prior = sdme)",
     {f0 <- frm(f_me, data = d)
      u <- update(f0, prior = set_prior("exponential(1)", class = "sdme"))
      as.data.frame(prior_summary(u))[, c("prior", "class")]})

cat("\n#### 7. cor group rows on the sample route (lane's convention claim)\n")
print(as.data.frame(default_prior(f_mere, data = d, route = "sample"))[
  , c("prior", "class", "coef", "group", "source")])
cat("DONE\n")
