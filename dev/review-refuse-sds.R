# Item 5: class sds on every path, on forms the lane did not use. Seed 9900.
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
set.seed(9900)
n <- 200
d <- data.frame(x = runif(n), z = runif(n), u = runif(n), w = rnorm(n),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(sample(1:8, n, TRUE)))
d$y <- sin(6 * d$x) + d$z^2 + rnorm(8, 0, 0.4)[d$g] + rnorm(n, 0, 0.4)
d$y2 <- cos(6 * d$u) + rnorm(n, 0, 0.5)
d$yp <- rpois(n, exp(0.5 + sin(3 * d$x)))
key <- function(t) {
  t <- as.data.frame(t)
  paste(t$class, t$coef, t$group, t$resp, t$dpar, t$nlpar, sep = "|")
}
tabcmp <- function(lab, form, fam = NULL) {
  fa <- eval(parse(text = form), envir = asNamespace("frmtmb"))
  fb <- eval(parse(text = form), envir = asNamespace("brms"))
  a <- tr(if (is.null(fam)) default_prior(fa, data = d) else default_prior(fa, data = d, family = eval(parse(text = fam), envir = asNamespace("frmtmb"))))
  b <- tr(if (is.null(fam)) brms::default_prior(fb, data = d) else brms::default_prior(fb, data = d, family = eval(parse(text = fam), envir = asNamespace("brms"))))
  if (is_err(a) || is_err(b)) { pr("TAB %-30s frm %s | brms %s", lab, msg(a, 90), msg(b, 90)); return(invisible()) }
  ka <- key(a); kb <- key(b)
  sa <- ka[grepl("^sds", ka)]; sb <- kb[grepl("^sds", kb)]
  pr("TAB %-30s set %s order %s sds-equal %s; frm-only [%s] brms-only [%s]", lab, setequal(ka, kb), identical(ka, kb), identical(sa, sb),
     paste(setdiff(ka, kb), collapse = " ; "), paste(setdiff(kb, ka), collapse = " ; "))
}
tabcmp("s(x, bs='cc')", "bf(y ~ s(x, bs = 'cc', k = 6))")
tabcmp("te(x,z)", "bf(y ~ te(x, z, k = 4))")
tabcmp("ti(x)+ti(z)+ti(x,z)", "bf(y ~ ti(x, k = 4) + ti(z, k = 4) + ti(x, z, k = 3))")
tabcmp("t2 by", "bf(y ~ f + t2(x, z, by = f, k = 4))")
tabcmp("s x z 2d", "bf(y ~ s(x, z, k = 8))")
tabcmp("s by g", "bf(y ~ g + s(x, by = g, k = 4))")
tabcmp("s bs=re", "bf(y ~ s(g, bs = 're'))")
tabcmp("s long deparse", "bf(y ~ s(x, bs = 'cr', k = 5, by = f, id = 1, m = 2) + f)")
tabcmp("mu+sigma", "bf(y ~ s(x, k = 5), sigma ~ s(z, k = 5))")
tabcmp("nl", "bf(y ~ a + b, a ~ s(x, k = 5), b ~ 1 + (1 | g), nl = TRUE)")
tabcmp("mv", "mvbf(bf(y ~ s(x, k = 5)), bf(y2 ~ s(u, k = 5)), rescor = FALSE)")
tabcmp("poly", "bf(y ~ s(I(x^2), k = 5))")
tabcmp("poisson", "bf(yp ~ s(x, k = 5))", "poisson()")
tabcmp("mix", "bf(y ~ s(x, k = 5))", "mixture(gaussian(), gaussian(), order = 'none')")
tabcmp("categorical", "bf(y3 ~ s(x, k = 5))", NULL) -> z

# validate_prior spellings and refusals, side by side
vp <- function(lab, form, prior, fam = NULL) {
  fa <- eval(parse(text = form), envir = asNamespace("frmtmb"))
  fb <- eval(parse(text = form), envir = asNamespace("brms"))
  a <- tr(validate_prior(eval(parse(text = prior), envir = asNamespace("frmtmb")), fa, data = d))
  b <- tr(brms::validate_prior(eval(parse(text = prior), envir = asNamespace("brms")), fb, data = d))
  pr("VP %-34s frm %-4s | brms %-4s | %s || %s", lab, if (is_err(a)) "ERR" else "ok", if (is_err(b)) "ERR" else "ok",
     if (is_err(a)) msg(a, 90) else "", if (is_err(b)) msg(b, 90) else "")
}
sp <- "set_prior('normal(0, 1)', class = 'sds'"
vp("class only, smooth", "bf(y ~ s(x))", paste0(sp, ")"))
vp("class only, no smooth", "bf(y ~ x)", paste0(sp, ")"))
vp("coef ok", "bf(y ~ s(x))", paste0(sp, ", coef = 's(x)')"))
vp("coef with spaces lost", "bf(y ~ s(x, k = 5))", paste0(sp, ", coef = 's(x,k=5)')"))
vp("coef spaced", "bf(y ~ s(x, k = 5))", paste0(sp, ", coef = 's(x, k = 5)')"))
vp("coef sx_1", "bf(y ~ s(x))", paste0(sp, ", coef = 'sx_1')"))
vp("coef other term", "bf(y ~ s(x))", paste0(sp, ", coef = 's(z)')"))
vp("group", "bf(y ~ s(x))", paste0(sp, ", group = 's(x)')"))
vp("group g", "bf(y ~ s(x) + (1 | g))", paste0(sp, ", group = 'g')"))
vp("bounds", "bf(y ~ s(x))", paste0(sp, ", lb = 0)"))
vp("ub", "bf(y ~ s(x))", paste0(sp, ", ub = 5)"))
vp("coef + lb", "bf(y ~ s(x))", paste0(sp, ", coef = 's(x)', lb = 0)"))
vp("dpar", "bf(y ~ s(x), sigma ~ s(z))", paste0(sp, ", dpar = 'sigma')"))
vp("dpar wrong", "bf(y ~ s(x), sigma ~ z)", paste0(sp, ", dpar = 'sigma')"))
vp("nlpar", "bf(y ~ a, a ~ s(x), nl = TRUE)", paste0(sp, ", nlpar = 'a')"))
vp("nlpar coef", "bf(y ~ a, a ~ s(x), nl = TRUE)", paste0(sp, ", nlpar = 'a', coef = 's(x)')"))
vp("resp", "mvbf(bf(y ~ s(x)), bf(y2 ~ s(u)), rescor = FALSE)", paste0(sp, ", resp = 'y')"))
vp("no resp, mv", "mvbf(bf(y ~ s(x)), bf(y2 ~ s(u)), rescor = FALSE)", paste0(sp, ")"))
vp("by-level coef", "bf(y ~ f + s(x, by = f))", paste0(sp, ", coef = 's(x, by = f)')"))
vp("by-level coef a", "bf(y ~ f + s(x, by = f))", paste0(sp, ", coef = 's(x,by=f):fa')"))
vp("t2 coef", "bf(y ~ t2(x, z))", paste0(sp, ", coef = 't2(x, z)')"))
vp("t2 coef 2", "bf(y ~ t2(x, z))", paste0(sp, ", coef = 't2xz_1')"))
sd <- "set_prior('normal(0, 1)', class = 'sd'"
vp("sd, smooth only", "bf(y ~ s(x))", paste0(sd, ")"))
vp("sd group smooth", "bf(y ~ s(x))", paste0(sd, ", group = 's(x)')"))
vp("sd group smooth + RE", "bf(y ~ s(x) + (1 | g))", paste0(sd, ", group = 's(x)')"))
vp("sd class, smooth + RE", "bf(y ~ s(x) + (1 | g))", paste0(sd, ")"))
vp("sd coef Intercept group g", "bf(y ~ s(x) + (1 | g))", paste0(sd, ", group = 'g', coef = 'Intercept')"))
vp("sd coef s(x)", "bf(y ~ s(x) + (1 | g))", paste0(sd, ", coef = 's(x)')"))
vp("sd no smooth no RE", "bf(y ~ x)", paste0(sd, ")"))
vp("sd dpar sigma, smooth", "bf(y ~ s(x), sigma ~ s(z))", paste0(sd, ", dpar = 'sigma')"))
vp("sd over nlpar smooth", "bf(y ~ a, a ~ s(x), nl = TRUE)", paste0(sd, ", nlpar = 'a')"))

# frm() and frm_sample() reach only the smooth: exact densities, by the taped
# value: with a prior on sds alone the penalized objective moves by exactly
# the prior at the fit's sd; a class sd prior beside it leaves sds alone
f0 <- frm(bf(y ~ s(x) + (1 | g)), data = d)
b0 <- f0$opt$par
gapf <- function(pri) {
  f <- frm(bf(y ~ s(x) + (1 | g)), data = d, prior = pri, start = NULL)
  f
}
f1 <- gapf(set_prior("exponential(2)", class = "sds"))
f2 <- gapf(set_prior("exponential(2)", class = "sd"))
pl <- function(f) {
  vc <- VarCorr(f)
  c(sm = unname(vc[["s(x)"]]$sd[1]), g = unname(vc[["g"]]$sd[1]))
}
pr("no prior sd %s; sds prior %s; sd prior %s", paste(round(pl(f0), 4), collapse = "/"),
   paste(round(pl(f1), 4), collapse = "/"), paste(round(pl(f2), 4), collapse = "/"))

# prior_summary, update, hypothesis, prior_draws, sample_prior on both routes
ps <- tr(prior_summary(f1))
print(as.data.frame(ps)[, c("prior", "class", "coef", "group")])
up <- tr(update(f1, . ~ . + w))
pr("update keeps the sds prior: %s; %s", msg(up, 100), if (is_err(up)) "" else {
   p2 <- as.data.frame(prior_summary(up)); paste(p2$class[nzchar(p2$prior)], collapse = ",")
})
up2 <- tr(update(f1, . ~ . - s(x)))
pr("update drops the smooth: %s", msg(up2, 150))
hy <- tr(hypothesis(f1, "sds_sx_1 > 0", class = NULL))
pr("hypothesis sds_sx_1: %s", msg(hy, 150))
hy2 <- tr(hypothesis(f1, "sd_g__Intercept > 0", class = NULL))
pr("hypothesis sd_g__Intercept: %s", msg(hy2, 150))
vn <- tr(variables(f1))
pr("variables: %s", if (is_err(vn)) msg(vn, 150) else paste(grep("sd", vn, value = TRUE), collapse = " "))
sm <- tr(frm_sample(bf(y ~ s(x, k = 5) + (1 | g)), data = d, prior = set_prior("exponential(2)", class = "sds"),
                    chains = 1, iter = 300, seed = 2, refresh = 0, sample_prior = "yes"))
pr("frm_sample sds + sample_prior: %s", msg(sm, 200))
if (!is_err(sm)) {
  pd <- tr(prior_draws(sm))
  pr("  prior_draws: %s; %s", msg(pd, 120), if (is_err(pd)) "" else paste(names(pd), collapse = " "))
  if (!is_err(pd)) {
    cn <- grep("sds", names(pd), value = TRUE)
    if (length(cn)) pr("  prior mean of %s %.3f (want 0.5)", cn[1], mean(pd[[cn[1]]]))
  }
  pr("  prior_summary rows: %s", paste(as.data.frame(prior_summary(sm))$class, collapse = ","))
  v <- tr(names(posterior::as_draws_df(sm)))
  pr("  draws: %s", paste(grep("^sd", v, value = TRUE), collapse = " "))
  h <- tr(hypothesis(sm, "sds_sx_1 > 0", class = NULL))
  pr("  hypothesis on draws: %s", msg(h, 140))
  h <- tr(hypothesis(sm, "sds_sx_1 > 0.1", class = NULL))
  pr("  hypothesis with prior draws: %s %s", msg(h, 140), if (is_err(h)) "" else format(h$hypothesis$Evid.Ratio, digits = 3))
}
# defaults on the sampling route against brms's stancode
fb <- brms::bf(y ~ s(x, k = 5) + (1 | g))
pb <- brms::default_prior(fb, data = d)
pa <- tr(default_prior(bf(y ~ s(x, k = 5) + (1 | g)), data = d, route = "frm_sample"))
pr("default_prior(route=): %s", msg(pa, 120))
print(as.data.frame(pb)[pb$class %in% c("sds", "sd"), c("prior", "class", "coef", "group", "lb")])
if (!is_err(pa)) print(as.data.frame(pa)[pa$class %in% c("sds", "sd"), c("prior", "class", "coef", "group", "lb")])

# frm_simulate_new / simulate-new reads priors too
cand <- ls(asNamespace("frmtmb"), pattern = "simulate_new|frm_simulate")
pr("simulate-new functions: %s", paste(cand, collapse = " "))
a <- tr(simulate_new(bf(y ~ s(x, k = 5)), newdata = d, prior = set_prior("exponential(2)", class = "sds"), seed = 3))
pr("simulate_new sds: %s", msg(a, 200))
a <- tr(frm_sbc <- NULL)

# brms's table through as_priorlist, with dpar and t2, and a by-smooth
for (fm in c("bf(y ~ t2(x, z, k = 4))", "bf(y ~ f + s(x, by = f, k = 4))", "bf(y ~ s(x, k = 5), sigma ~ s(z, k = 5))")) {
  tb <- brms::default_prior(eval(parse(text = fm), envir = asNamespace("brms")), data = d)
  tb$prior[tb$class == "sds" & !nzchar(tb$coef)] <- "normal(0, 1)"
  a <- tr(frm(eval(parse(text = fm), envir = asNamespace("frmtmb")), data = d, prior = tb))
  pr("brms table into frm(): %-44s %s", fm, msg(a, 140))
}

# smooth in dpar and nlpar reach, the theta class coefficient table and the
# one-slot refusal on the fit route
f3 <- tr(frm(bf(y ~ s(x, k = 5), sigma ~ s(z, k = 5)), data = d,
             prior = set_prior("exponential(2)", class = "sds", dpar = "sigma")))
f4 <- tr(frm(bf(y ~ s(x, k = 5), sigma ~ s(z, k = 5)), data = d))
if (!is_err(f3) && !is_err(f4)) {
  pr("dpar sds: sd mu %s vs %s, sigma %s vs %s", format(VarCorr(f3)[[1]]$sd[1], digits = 4), format(VarCorr(f4)[[1]]$sd[1], digits = 4),
     format(VarCorr(f3)[[2]]$sd[1], digits = 4), format(VarCorr(f4)[[2]]$sd[1], digits = 4))
} else pr("dpar sds: %s %s", msg(f3, 120), msg(f4, 120))
