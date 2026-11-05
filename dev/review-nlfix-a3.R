# Area 3: multi-mu spellings, brms stancode against frm(). Seed 5301.
source("dev/review-nlfix-common.R")
set.seed(5301)
n <- 200
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3), ntr = 10L)
d$cat <- factor(sapply(seq_len(n), function(i)
  sample(c("a", "b", "c"), 1, prob = c(1, exp(0.8 * d$x[i]), exp(-0.5 * d$x[i])))))
m <- t(vapply(seq_len(n), function(i)
  stats::rmultinom(1, 10, c(1, exp(0.8 * d$x[i]), exp(-0.5 * d$x[i])))[, 1], numeric(3)))
colnames(m) <- c("a", "b", "c")
d$cnt <- m
s <- (m + 0.5) / 11.5
d$sim <- s
d$yg <- ifelse(runif(n) < 0.4, -1, 1) + 0.5 * d$x + rnorm(n, 0, 0.5)
brm_code <- function(txt, fam) {
  r <- tryCatch(suppressWarnings(suppressMessages({
    f <- eval(parse(text = txt), envir = asNamespace("brms"))
    brms::stancode(f, data = d, family = eval(parse(text = fam), envir = asNamespace("brms")))
  })), error = function(e) e)
  if (inherits(r, "error")) return(paste("REFUSES:", gsub("\n", " ", conditionMessage(r))))
  ln <- strsplit(as.character(r), "\n")[[1]]
  keep <- grep("nlp_[a-zA-Z0-9_]*\\[n\\]|mu[a-zA-Z0-9_]*\\[n\\] =|mu[a-zA-Z0-9_]* =|mu[a-zA-Z0-9_]* \\+=|target \\+= [a-z_]+_l", ln, value = TRUE)
  paste(unique(trimws(keep)), collapse = " || ")
}
frm_try <- function(txt, fam, ...) {
  r <- tt(frm(eval(parse(text = txt)), family = eval(parse(text = fam)), data = d, ...))
  if (iserr(r)) return(paste("REFUSES:", substr(r$msg, 1, 300)))
  paste0("ll=", format(as.numeric(logLik(r)), digits = 9), " beta=",
         paste(names(r$estimates$beta), collapse = ","), " betad=",
         paste(names(r$estimates$betad), collapse = ","),
         if (length(attr(r, "w"))) paste0(" WARN[", substr(attr(r, "w")[1], 1, 100), "]"))
}
specs <- list(
  shared_cat = c("bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", "categorical()"),
  perloc_cat = c("bf(cat ~ 1, nlf(mub ~ ab + bb * x), nlf(muc ~ ac + bc * x), ab + bb + ac + bc ~ 1)", "categorical()"),
  perloc_cat2 = c("bf(cat ~ 1) + nlf(mub ~ ab + bb * x) + nlf(muc ~ ac + bc * x) + lf(ab + bb + ac + bc ~ 1)", "categorical()"),
  perloc_flag = c("bf(cat ~ ab + bb * x, muc ~ x, ab + bb ~ 1, nl = TRUE)", "categorical()"),
  one_nlf_lin = c("bf(cat ~ x, nlf(muc ~ ac + bc * x), ac + bc ~ 1)", "categorical()"),
  one_nlf_flag = c("bf(cat ~ 1, nlf(mub ~ ab + bb * x), ab + bb ~ 1, nl = TRUE)", "categorical()"),
  mu2_b = c("bf(cat ~ 1, nlf(mub ~ mub_a + mub_b * x), nlf(muc ~ muc_a + muc_b * x), mub_a + mub_b + muc_a + muc_b ~ 1)", "categorical()"),
  perloc_mn = c("bf(cnt | trials(ntr) ~ 1, nlf(mub ~ ab + bb * x), nlf(muc ~ ac + bc * x), ab + bb + ac + bc ~ 1)", "multinomial()"),
  perloc_dir = c("bf(sim ~ 1, nlf(muV2 ~ ab + bb * x), nlf(muV3 ~ ac + bc * x), ab + bb + ac + bc ~ 1)", "dirichlet()"),
  perloc_ln = c("bf(sim ~ 1, nlf(muV2 ~ ab + bb * x), nlf(muV3 ~ ac + bc * x), ab + bb + ac + bc ~ 1)", "logistic_normal()"),
  shared_mix = c("bf(yg ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", "mixture(gaussian(), gaussian())"),
  shared_mix_ord = c("bf(yg ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", "mixture(gaussian(), gaussian(), order = 'mu')"),
  perloc_mix = c("bf(yg ~ 1, nlf(mu1 ~ a1 + b1 * x), nlf(mu2 ~ a2 + b2 * x), a1 + b1 + a2 + b2 ~ 1)", "mixture(gaussian(), gaussian())"),
  perloc_mix_ord = c("bf(yg ~ 1, nlf(mu1 ~ a1 + b1 * x), nlf(mu2 ~ a2 + b2 * x), a1 + b1 + a2 + b2 ~ 1)", "mixture(gaussian(), gaussian(), order = 'mu')"),
  perloc_mix_mu2_b = c("bf(yg ~ 1, nlf(mu1 ~ mu1_a + mu1_b * x), nlf(mu2 ~ mu2_a + mu2_b * x), mu1_a + mu1_b + mu2_a + mu2_b ~ 1)", "mixture(gaussian(), gaussian())"),
  mix_one_lin = c("bf(yg ~ x, nlf(mu2 ~ a2 + b2 * x), a2 + b2 ~ 1)", "mixture(gaussian(), gaussian())"),
  mix_flag_lin = c("bf(yg ~ a + b * x, mu2 ~ x, a + b ~ 1, nl = TRUE)", "mixture(gaussian(), gaussian())"),
  mix_theta = c("bf(yg ~ b0 + b1 * x, theta1 ~ z, b0 + b1 ~ 1, nl = TRUE)", "mixture(gaussian(), gaussian())"),
  mix_sigma = c("bf(yg ~ b0 + b1 * x, sigma2 ~ z, b0 + b1 ~ 1, nl = TRUE)", "mixture(gaussian(), gaussian())"),
  mix_pois = c("bf(ntr ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", "mixture(poisson(), poisson())"),
  mix_cum = c("bf(ord ~ b1 * x, b1 ~ 1, nl = TRUE)", "mixture(cumulative(), cumulative())")
)
d$ord <- cut(d$yg, c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
d$ntr <- rpois(n, exp(1 + 0.3 * d$x))
for (nm in names(specs)) {
  sp <- specs[[nm]]
  cat("\n==", nm, ":", sp[1], "|", sp[2], "\n")
  cat("  brms:", substr(brm_code(sp[1], sp[2]), 1, 700), "\n")
  cat("  frm :", substr(frm_try(sp[1], sp[2]), 1, 500), "\n")
}
cat("\n=== priors per nlpar under multi-mu\n")
f <- bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)
pf <- as.data.frame(tt(get_prior(f, data = d, family = categorical())))
pb <- as.data.frame(tt(brms::default_prior(brms::bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), data = d, family = brms::categorical())))
print(pf[, intersect(c("prior", "class", "coef", "group", "dpar", "nlpar"), names(pf))])
print(pb[, intersect(c("prior", "class", "coef", "group", "dpar", "nlpar"), names(pb))])
cat("MAP with a prior on nlpar b1 under categorical, tight at 0 vs without\n")
r0 <- tt(frm(f, family = categorical(), data = d))
r1 <- tt(frm(f, family = categorical(), data = d, prior = set_prior("normal(0, 0.01)", nlpar = "b1")))
show("r0", if (iserr(r0)) r0 else r0$estimates$beta); show("r1", if (iserr(r1)) r1 else r1$estimates$beta)
r2 <- tt(frm(f, family = categorical(), data = d, prior = set_prior("normal(0, 0.01)", nlpar = "b1", dpar = "mub")))
show("r2 (dpar = mub)", if (iserr(r2)) r2 else r2$estimates$beta)
fb <- tryCatch(brms::validate_prior(brms::set_prior("normal(0, 0.01)", nlpar = "b1", dpar = "mub"),
     brms::bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), data = d, family = brms::categorical()), error = function(e) conditionMessage(e))
cat("brms dpar= mub on an nlpar:", if (is.character(fb)) substr(fb, 1, 300) else "accepted", "\n")
cat("\n=== refcat under nl\n")
for (rc in c("'a'", "'c'", "NA")) {
  fam <- paste0("categorical(refcat = ", rc, ")")
  cat(fam, "\n  brms:", substr(brm_code("bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", fam), 1, 400), "\n")
  cat("  frm :", frm_try("bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", fam), "\n")
}
cat("\n=== readers on a shared categorical body\n")
fit <- r0
for (nm in c("fitted", "predict")) {
  x <- tt(get(nm)(fit, newdata = d[1:5, ])); show(nm, if (iserr(x)) x else dim(x))
}
x <- tt(frm_linpred(fit, dpar = "mub")); show("linpred mub", if (iserr(x)) x else range(x))
x <- tt(frm_linpred(fit, dpar = "muc")); show("linpred muc", if (iserr(x)) x else range(x))
x <- tt(frm_curve(fit, newdata = d[1:5, ])); show("frm_curve", if (iserr(x)) x else nrow(x))
x <- tt(simulate(fit, nsim = 2)); show("simulate", if (iserr(x)) x else length(x))
x <- tt(confint(fit)); show("confint", if (iserr(x)) x else nrow(x))
x <- tt(update(fit, data = d[1:150, ])); show("update", if (iserr(x)) x else as.numeric(logLik(x)))
x <- tt(conditional_effects(fit, effects = "x", categorical = TRUE)); show("cond eff", if (iserr(x)) x else 1)
x <- tt(marginaleffects::predictions(fit, newdata = d[1:5, ])); show("predictions", if (iserr(x)) x else nrow(x))
x <- tt(summary(fit)); print(x)
