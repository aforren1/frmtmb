# Review of lane meme, short check of fix round 1: frmtmb side.
#   Rscript dev/review-meme-sc-frm.R > dev/review-meme-out/sc/frm.txt
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-meme"
source(file.path(wt, "dev/review-meme-shapes.R"))
source(file.path(wt, "dev/review-meme-sc-shapes.R"))
d <- sc_data()
sh <- sc_shapes(bf, mvbf, set_rescor, cumulative, gp)
B <- readRDS(file.path(wt, "dev/review-meme-out/sc/brms.rds"))
q <- function(expr) suppressMessages(suppressWarnings(expr))
try_ <- function(lab, expr) {
  r <- tryCatch(q(expr), error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("--", lab, "\n")
  if (is.character(r) && length(r) == 1) cat("  ", r, "\n") else print(r)
  invisible(r)
}
key <- function(t) {
  t <- as.data.frame(t)
  lb <- suppressWarnings(as.numeric(as.character(t$lb)))
  ub <- suppressWarnings(as.numeric(as.character(t$ub)))
  paste(ifelse(nzchar(t$prior), t$prior, "(flat)"), t$class, t$coef,
        t$group, t$resp, t$dpar, ifelse(is.na(lb), "", lb),
        ifelse(is.na(ub), "", ub), sep = "|")
}
ents <- function(fit, pr) {
  rs <- q(frmtmb.sample:::sample_resolve_priors(fit, pr))
  vapply(rs$ri$entries, function(e) {
    paste0(e$comp, "[", paste(e$idx, collapse = ","), "]:", e$dist$kind,
           e$dist$eta %||% e$dist$rate %||% e$dist$scale %||% "")
  }, "")
}

cat("#### A. sampling-route tables against brms's live default_prior()\n")
nsame <- 0
for (nm in names(sh)) {
  s <- sh[[nm]]
  ft <- tryCatch(default_prior(s$f, data = d, family = s$fam,
                               route = "sample"),
                 error = function(e) conditionMessage(e))
  b <- B[[nm]]
  if (is.character(ft) || is.character(b)) {
    cat("==", nm, "frmtmb:", if (is.character(ft)) ft else "table",
        " brms:", if (is.character(b)) b else "table", "\n")
    next
  }
  kb <- key(b)
  kf <- key(ft)
  same <- identical(kb, kf)
  nsame <- nsame + same
  cat("==", nm, if (same) "IDENTICAL" else "DIFFERENT", "\n")
  if (!same) {
    ob <- setdiff(kb, kf)
    of <- setdiff(kf, kb)
    if (length(ob)) cat("   brms only:  ", ob, sep = "\n     ")
    if (length(of)) cat("   frmtmb only:", of, sep = "\n     ")
    if (!length(ob) && !length(of)) cat("   order differs\n")
    cat("\n")
  }
}
cat(sprintf("RESULT tables identical %d of %d\n", nsame, length(sh)))

cat("\n#### B. group-row round trips on the sampling route\n")
for (nm in c("thres", "thres_eq")) {
  s <- sh[[nm]]
  fit0 <- q(frm(s$f, data = d, family = s$fam, dry_run = "objective"))
  ts <- default_prior(s$f, data = d, family = s$fam, route = "sample")
  print(as.data.frame(ts)[ts$class %in% c("Intercept", "delta"),
                          c("prior", "class", "coef", "group", "lb")])
  cat("  as printed:", ents(fit0, ts), "\n", fill = 78)
  t1 <- ts
  cls <- if (nm == "thres") "Intercept" else "delta"
  t1$prior[t1$class == cls & !nzchar(t1$group) & !nzchar(t1$coef)] <-
    if (nm == "thres") "normal(0, 7)" else "exponential(3)"
  cat("  class row edited:", ents(fit0, t1), "\n", fill = 78)
  t2 <- t1
  t2$prior[t2$class == cls & t2$group == "b" & !nzchar(t2$coef)] <-
    if (nm == "thres") "normal(0, 9)" else "exponential(9)"
  cat("  class row and group b edited:", ents(fit0, t2), "\n", fill = 78)
}
s <- sh$reme
fit0 <- q(frm(s$f, data = d, dry_run = "objective"))
ts <- default_prior(s$f, data = d, route = "sample")
t3 <- ts
t3$prior[t3$class == "corme" & !nzchar(t3$group)] <- "lkj(2)"
t3$prior[t3$class == "corme" & t3$group == "g"] <- "lkj(5)"
t3$prior[t3$class == "cor" & !nzchar(t3$group)] <- "lkj(3)"
t3$prior[t3$class == "sd" & !nzchar(t3$group) & !nzchar(t3$coef)] <-
  "exponential(3)"
t3$prior[t3$class == "sd" & t3$group == "g" & !nzchar(t3$coef)] <-
  "exponential(7)"
cat("  reme, class rows edited and explicit group rows:", ents(fit0, t3),
    "\n", fill = 78)
cat("  reme, set_prior class + group:",
    ents(fit0, set_prior("lkj(2)", class = "corme") +
           set_prior("lkj(5)", class = "corme", group = "g") +
           set_prior("exponential(3)", class = "sd") +
           set_prior("exponential(7)", class = "sd", group = "g")),
    "\n", fill = 78)

cat("\n#### C. frm()'s route\n")
tf <- default_prior(s$f, data = d)
print(as.data.frame(tf)[tf$class %in% c("sd", "cor", "corme"),
                        c("prior", "class", "coef", "group")])
try_("frm(prior = sample table) equals hand-written", {
  a <- frm(s$f, data = d, prior = ts)
  b <- frm(s$f, data = d, prior =
             set_prior("student_t(3, 0, 2.5)", class = "sd") +
             set_prior("lkj(1)", class = "cor") +
             set_prior("lkj(1)", class = "corme") +
             set_prior(ts$prior[ts$class == "Intercept"], class = "Intercept") +
             set_prior("student_t(3, 0, 2.5)", class = "sigma"))
  c(table = as.numeric(logLik(a)), hand = as.numeric(logLik(b)))
})

cat("\n#### D. coef-bound refusal on frm_sample()\n")
f <- bf(y ~ x + w + (1 + x | g) + me(xg, sg, gr = g) + me(wg, sg, gr = g))
fitf <- q(frm(f, data = d, dry_run = "objective"))
ps <- list(
  b_coef_ub = set_prior("normal(0, 1)", class = "b", coef = "w", ub = 0.1),
  sd_coef_ub = set_prior("normal(0, 1)", class = "sd", group = "g",
                         coef = "x", ub = 2),
  sd_coef_lb0 = set_prior("normal(0, 1)", class = "sd", group = "g",
                          coef = "x", lb = 0),
  sd_coef_lb0_ub2 = set_prior("normal(0, 1)", class = "sd", group = "g",
                              coef = "x", lb = 0, ub = 2),
  sd_coef_lb05 = set_prior("normal(0, 1)", class = "sd", group = "g",
                           coef = "x", lb = 0.5),
  sd_group_ub = set_prior("normal(0, 1)", class = "sd", group = "g", ub = 2),
  meanme_coef_lb = set_prior("normal(0, 1)", class = "meanme", coef = "mexg",
                             lb = 0),
  sdme_coef_ub = set_prior("exponential(1)", class = "sdme", coef = "mexg",
                           ub = 0.5),
  sdme_coef_lb0 = set_prior("exponential(1)", class = "sdme", coef = "mexg",
                            lb = 0),
  sdme_class_ub = set_prior("exponential(1)", class = "sdme", ub = 3),
  b_dpar_coef_lb = set_prior("normal(0, 1)", class = "b", coef = "x",
                             dpar = "sigma", lb = 0))
for (k in names(ps)) {
  ff <- if (k == "b_dpar_coef_lb") {
    q(frm(bf(y ~ x, sigma ~ x), data = d, dry_run = "objective"))
  } else fitf
  r <- tryCatch({ents(ff, ps[[k]]); "accepted"},
                error = function(e) paste("ERROR:", substr(conditionMessage(e),
                                                           1, 70)))
  cat(sprintf("%-16s %s\n", k, r))
}
cat("-- frm() keeps a coef bound:\n")
try_("frm b coef w ub 0.1", {
  ft <- frm(f, data = d, prior = ps$b_coef_ub)
  ft$estimates$beta[["w"]] %||% fixef(ft)["w", 1]
})
try_("frm_sample(<MAP fit with a coef bound>)", {
  ft <- frm(f, data = d, prior = ps$b_coef_ub)
  frm_sample(ft, chains = 1, iter = 50, refresh = 0)
  "sampled"
})
cat("-- exemptions\n")
fn <- bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE)
dn <- data.frame(x = runif(60, 0, 3))
set.seed(2)
dn$y <- 2 * exp(-0.7 * dn$x) + rnorm(60, 0, 0.1)
pn <- set_prior("normal(1, 1)", class = "Intercept", nlpar = "a", lb = 0) +
  set_prior("normal(1, 1)", class = "Intercept", nlpar = "b", lb = 0.4)
try_("nlpar Intercept with lb, prior draws", {
  ds <- frm_sample(fn, data = dn, family = gaussian(), prior = pn +
                     set_prior("exponential(1)", class = "sigma"),
                   sample_prior = "yes", chains = 2, iter = 600, seed = 3,
                   refresh = 0)
  pd <- prior_draws(ds)
  m <- as.matrix(ds)
  list(names = names(pd),
       prior_min = vapply(pd[grep("a|b", names(pd))], min, 0),
       post_min = apply(m[, grep("^b_[ab]_", colnames(m)), drop = FALSE], 2,
                        min))
})
try_("brms-spelled nlpar intercept: b coef Intercept nlpar b lb 0.4",
     ents(q(frm(fn, data = dn, dry_run = "objective")),
          set_prior("normal(1, 1)", class = "b", coef = "Intercept",
                    nlpar = "b", lb = 0.4)))
try_("theta coef with a bound (ar)", {
  fa <- bf(y ~ x + ar(time = tt, gr = g, p = 1))
  da <- d
  da$tt <- rep(1:5, 8)
  fa0 <- frm(fa, data = da, dry_run = "objective")
  dp <- as.data.frame(default_prior(fa, data = da))
  print(dp[dp$class == "theta" | dp$class == "ar", c("prior", "class", "coef", "lb", "ub")])
  ents(fa0, set_prior("", class = "ar", lb = -0.5, ub = 0.5))
})
try_("gp table passed back to frm_sample", {
  tg <- default_prior(sh$gp$f, data = d, route = "sample")
  print(as.data.frame(tg)[, c("prior", "class", "coef", "lb", "ub")])
  ds <- frm_sample(sh$gp$f, data = d, family = gaussian(), prior = tg,
                   chains = 1, iter = 100, seed = 1, refresh = 0)
  "sampled"
})
try_("gp table, coef row lb edited to 0.1", {
  tg <- default_prior(sh$gp$f, data = d, route = "sample")
  tg$lb[tg$class == "sdgp" & nzchar(tg$coef)] <- 0.1
  frm_sample(sh$gp$f, data = d, family = gaussian(), prior = tg,
             chains = 1, iter = 100, seed = 1, refresh = 0)
  "sampled"
})

cat("\n#### E. me(x) without sdx\n")
fs <- bf(y ~ me(x) + z)
ft <- default_prior(fs, data = d)
cat("fit route vs brms keys (no prior column):",
    identical(sub("^[^|]*\\|", "", key(ft)), sub("^[^|]*\\|", "", key(B$nosdx))),
    "\n")
vp <- validate_prior(set_prior("normal(0, 1)", class = "meanme"), fs, data = d)
print(as.data.frame(vp)[, c("prior", "class", "coef", "lb", "source")])
try_("frm()", frm(fs, data = d))
try_("frm_sample()", frm_sample(fs, data = d, family = gaussian(), chains = 1,
                                iter = 50, refresh = 0))
try_("frm_sample(sample_prior = 'only')",
     frm_sample(fs, data = d, family = gaussian(), sample_prior = "only",
                prior = set_prior("normal(0,1)", class = "b"), chains = 1,
                iter = 50, refresh = 0))
cat("DONE\n")
