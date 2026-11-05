# Short check, items 4, 5, 6, 7. Seed 11300.
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
set.seed(11300)
quiet <- function(e) suppressMessages(suppressWarnings(e))
n <- 90; x <- rnorm(n)
d <- data.frame(x = x, y = rgamma(n, 2, 2 / exp(0.2 + 0.4 * x)), w = rnorm(n), z = rnorm(n), s = 0.3, Y_a = rgamma(n, 3, 2))
M <- seq(4, n, by = 15); M2 <- seq(7, n, by = 15)
d$y[M] <- NA; d$w[M2] <- NA; d$Y_a[M2] <- NA
d$z <- 0.3 + 0.5 * ifelse(is.na(d$y), 1, d$y) + rnorm(n)
Gam <- stats::Gamma(link = "log"); Gb <- brms::brmsfamily("gamma", "log")

# 5. the names
pars_of <- function(form, data, data2 = NULL) {
  pri <- E$brms_flat_prior(form, data = data, family = NULL)
  code <- brms::make_stancode(form, data = data, prior = pri)
  sdat <- E$brms_standata(form, data = data, prior = pri)
  sf <- quiet(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 1, iter = 60, refresh = 0, seed = 1))
  grep("^(Ymi|Yl)", sf@sim$fnames_oi, value = TRUE)
}
cmp_names <- function(lab, ff, fb, data) {
  s <- tr(quiet(frm_sample(ff, data = data, chains = 1, iter = 60, seed = 1, refresh = 0)))
  if (is_err(s)) { pr("%-28s frm_sample: %s", lab, msg(s, 150)); return(invisible()) }
  a <- grep("^(Ymi|Yl)", names(posterior::as_draws_df(s)), value = TRUE)
  b <- tr(pars_of(fb, data))
  if (is_err(b)) { pr("%-28s brms: %s", lab, msg(b, 150)); return(invisible(s)) }
  pr("%-28s frm %d names, brms %d; same set %s; same order %s; frm head %s; brms head %s", lab, length(a), length(b), setequal(a, b), identical(a, b),
     paste(head(a, 2), collapse = " "), paste(head(b, 2), collapse = " "))
  if (!setequal(a, b)) pr("   only frm: %s; only brms: %s", paste(head(setdiff(a, b), 4), collapse = " "), paste(head(setdiff(b, a), 4), collapse = " "))
  invisible(s)
}
s1 <- cmp_names("univariate Ymi", bf(y | mi() ~ x) + Gam, brms::bf(y | mi() ~ x, family = Gb), d)
cmp_names("two mi responses", bf(y | mi() ~ x) + Gam + bf(w | mi() ~ x) + gaussian(),
          brms::bf(y | mi() ~ x, family = Gb) + brms::bf(w | mi() ~ x) + brms::set_rescor(FALSE), d)
cmp_names("mv with a reader", bf(y | mi() ~ x) + Gam + bf(z ~ mi(y) + x) + gaussian(),
          brms::bf(y | mi() ~ x, family = Gb) + brms::bf(z ~ mi(y) + x) + brms::set_rescor(FALSE), d)
cmp_names("resp with underscore", bf(Y_a | mi() ~ x) + Gam + bf(z ~ x) + gaussian(),
          brms::bf(Y_a | mi() ~ x, family = Gb) + brms::bf(z ~ x) + brms::set_rescor(FALSE), d)
cmp_names("mi(s) univariate", bf(y | mi(s) ~ x) + Gam, brms::bf(y | mi(s) ~ x, family = Gb), d)
cmp_names("mi(s) mv", bf(y | mi(s) ~ x) + Gam + bf(z ~ mi(y)) + gaussian(),
          brms::bf(y | mi(s) ~ x, family = Gb) + brms::bf(z ~ mi(y)) + brms::set_rescor(FALSE), d)
cmp_names("mi(s) and mi()", bf(y | mi(s) ~ x) + Gam + bf(w | mi() ~ x) + gaussian(),
          brms::bf(y | mi(s) ~ x, family = Gb) + brms::bf(w | mi() ~ x) + brms::set_rescor(FALSE), d)

# which row Ymi[k] is: value against the row's fitted value in brms's order,
# with the read and unread mixed on the fit route
dm <- d; dm$ok <- TRUE; dm$ok[M[1:3]] <- FALSE; dm$z2 <- ifelse(dm$ok, dm$z, 0); dm$id <- seq_len(n)
fr <- tr(frm(bf(y | mi() + index(id) ~ x) + Gam + bf(z2 | subset(ok) ~ mi(y, idx = id) + x) + gaussian(), data = dm))
pr("read/unread mix, fit: %s; latent %s", msg(fr, 150), if (is_err(fr)) "-" else length(fr$estimates$miss))
if (!is_err(fr)) {
  pr("  fit's labels: %s", paste(frmtmb:::brms_miss_labels(fr), collapse = " "))
  sf <- tr(quiet(frm_sample(fr, chains = 1, iter = 100, seed = 1, refresh = 0)))
  pr("  frm_sample(fit): %s", msg(sf, 150))
  if (!is_err(sf)) {
    dd <- posterior::as_draws_df(sf); nn <- grep("^Ymi", names(dd), value = TRUE)
    pr("  draws: %d of Ymi: %s", length(nn), paste(head(nn, 8), collapse = " "))
    # the first three are unread, the rest read: the read ones sit near the fit's
    pr("  value range %.3f %.3f (positive: Gamma)", min(as.matrix(as.data.frame(dd)[nn])), max(as.matrix(as.data.frame(dd)[nn])))
  }
}
# readers of the names
s <- s1
for (fn in c("variables", "as_draws_df", "posterior_summary", "summary", "fixef", "ranef", "prior_summary")) {
  v <- tr(do.call(switch(fn, as_draws_df = posterior::as_draws_df, posterior_summary = posterior_summary, get(fn)), list(s)))
  pr("  %-18s %s %s", fn, msg(v, 100), if (is_err(v)) "" else {
    nm <- if (fn %in% c("variables")) v else if (fn == "as_draws_df") names(v) else rownames(v)
    paste(head(grep("Ymi|miss", nm, value = TRUE), 2), collapse = " ")
  })
}
pp <- tr(quiet(prior_draws(s))); pr("  prior_draws without sample_prior: %s", msg(pp, 130))
up <- tr(quiet(update(s, . ~ . , chains = 1, iter = 60, refresh = 0))); pr("  update: %s", msg(up, 150))
pc <- tr(quiet(pp_check(s, ndraws = 3))); pr("  pp_check: %s", msg(pc, 130))
ce <- tr(quiet(conditional_effects(s))); pr("  conditional_effects: %s", msg(ce, 130))
h <- tr(hypothesis(s, "Ymi[1] > 0", class = NULL)); pr("  hypothesis on Ymi[1]: %s", msg(h, 130))
hb <- tr(posterior::subset_draws(posterior::as_draws_df(s), "Ymi[1]")); pr("  subset_draws Ymi[1]: %s", msg(hb, 100))
ss <- tr(quiet(frm_sample(bf(y | mi() ~ x) + Gam, data = d, chains = 1, iter = 100, seed = 1, refresh = 0, sample_prior = "yes")))
pr("  sample_prior yes: %s", msg(ss, 150))
if (!is_err(ss)) { pd <- tr(prior_draws(ss)); pr("  prior_draws: %s", if (is_err(pd)) msg(pd, 100) else paste(names(pd), collapse = " ")) }
pr("  mi_sampling exported: %s; help alias %s", "mi_sampling" %in% getNamespaceExports("frmtmb"), length(utils::help("mi_sampling", package = "frmtmb")))
a <- frm(bf(y | mi() ~ x) + Gam + bf(z ~ mi(y)) + gaussian(), data = d)
v <- tr(mi_sampling(a, "bad")); pr("  mi_sampling(what = 'bad'): %s", msg(v, 100))
v <- tr(mi_sampling(a, "natural")); pr("  natural: %s", if (is_err(v)) msg(v, 100) else paste(head(v$names, 3), collapse = " "))
b <- mi_sampling(a); pr("  fit form: latent %d; second call changes nothing: %s", length(b$estimates$miss), identical(length(mi_sampling(b)$estimates$miss), length(b$estimates$miss)))
pr("  logLik of the restored fit %.4f vs the fit %.4f; methods on it: fitted %s", as.numeric(logLik(b)), as.numeric(logLik(a)), msg(tr(fitted(b)), 80))

# 4. refusals, both orders; absent cases
set.seed(11301)
n <- 120; lev <- paste0("L", 1:8)
dm <- data.frame(g1 = sample(lev, n, TRUE), g2 = sample(lev, n, TRUE), x = rnorm(n), y = rnorm(n), w1 = runif(n), w2 = runif(n), v1 = runif(n), v2 = runif(n))
dm$W <- cbind(dm$w1, dm$w2)
fr <- function(lab, txt) {
  a <- tr(frm(eval(parse(text = txt), envir = asNamespace("frmtmb")), data = dm, dry_run = "frame"))
  b <- tr(brms::make_standata(eval(parse(text = txt), envir = asNamespace("brms")), data = dm))
  pr("%-38s frm %-60s brms %s", lab, msg(a, 58), msg(b, 40))
}
fr("key, weights differ", "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(v1, v2))))")
fr("key, weights differ, swapped", "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(v1, v2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2))))")
fr("key, weights then none", "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2)))")
fr("key, none then weights", "bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2))))")
fr("key, scale differs", "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2), scale = FALSE)), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2))))")
fr("key, same weights: fits", "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2))))")
fr("key, same, other spacing", "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1,w2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind( w1 , w2 ))))")
fr("key, same weights by a matrix", "bf(y ~ x + (1 | p | mm(g1, g2, weights = W)), sigma ~ 1 + (1 | p | mm(g1, g2, weights = W)))")
fr("key, W and cbind of the same", "bf(y ~ x + (1 | p | mm(g1, g2, weights = W)), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2))))")
fr("key, none and equal weights", "bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2)))")
fr("two keys, weights differ: fits", "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | q | mm(g1, g2, weights = cbind(v1, v2))))")
fr("no key, weights differ: fits", "bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | mm(g1, g2, weights = cbind(v1, v2))))")
dm$g <- factor(dm$g1)
fr("cor F + T", "bf(y ~ x + (1 | gr(g, id = 'p', cor = FALSE)), sigma ~ (1 | gr(g, id = 'p')))")
fr("cor T + F", "bf(y ~ x + (1 | gr(g, id = 'p')), sigma ~ (1 | gr(g, id = 'p', cor = FALSE)))")
fr("|p| spelled, T + F", "bf(y ~ x + (1 | p | g), sigma ~ (1 | p | gr(g, cor = FALSE)))")
fr("cor F + F: fits", "bf(y ~ x + (1 | gr(g, id = 'p', cor = FALSE)), sigma ~ (1 | gr(g, id = 'p', cor = FALSE)))")
fr("one term cor F: fits", "bf(y ~ x + (1 + x | gr(g, id = 'p', cor = FALSE)))")
fr("other key cor T: fits", "bf(y ~ x + (1 | gr(g, id = 'p', cor = FALSE)), sigma ~ (1 | gr(g, id = 'q')))")
fr("same label, other group: fits", "bf(y ~ x + (1 | gr(g, id = 'p', cor = FALSE)), sigma ~ (1 | gr(g1, id = 'p')))")
fr("mm cor F + T", "bf(y ~ x + (1 | mm(g1, g2, id = 'p', cor = FALSE)), sigma ~ (1 | mm(g1, g2, id = 'p')))")
fr("|| beside a key, mixed", "bf(y ~ x + (1 + x || g), sigma ~ (1 | p | g))")
# m2's tolerance
A0 <- exp(-abs(outer(1:8, 1:8, "-")) / 3); dimnames(A0) <- list(lev, lev)
tryA <- function(lab, A) {
  a <- tr(frm(bf(y ~ x + (1 | gr(g1, cov = A))), data = dm, data2 = list(A = A), dry_run = "frame"))
  b <- tr(brms::make_standata(brms::bf(y ~ x + (1 | gr(g1, cov = A))), data = dm, data2 = list(A = A)))
  pr("A %-34s frm %-50s brms %s", lab, msg(a, 48), msg(b, 48))
}
set.seed(1)
noise <- function(e, s = 1) { N <- matrix(rnorm(64, 0, e), 8, 8); A0 * s + N - t(N) * 0 }
for (e in c(1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-10, 1e-8)) {
  A <- A0; A[upper.tri(A)] <- A[upper.tri(A)] * (1 + rnorm(28, 0, e))
  tryA(sprintf("rel. rounding %g", e), A)
}
for (sc in c(1e-8, 1e8)) {
  A <- A0 * sc; A[1, 2] <- A[1, 2] * (1 + 1e-14)
  tryA(sprintf("scale %g, rounding 1e-14", sc), A)
  A <- A0 * sc; A[1, 2] <- A[1, 2] * 1.01
  tryA(sprintf("scale %g, a real 1%%", sc), A)
}
A <- crossprod(matrix(rnorm(64), 8, 8)); dimnames(A) <- list(lev, lev)
tryA("crossprod: exact", A)
K <- tcrossprod(matrix(rnorm(160), 8, 20)); K2 <- K %*% solve(K) %*% K; dimnames(K2) <- list(lev, lev)
pr("  asymmetry of a computed K %.1e", max(abs(K2 - t(K2))) / max(abs(K2)))
tryA("computed K %*% solve(K) %*% K", K2)
Ap <- A0; dimnames(Ap) <- list(lev, lev); Ax <- A0; Ax[1, 1] <- 1e-12
tryA("barely PD", Ax)
As <- A0; ev <- eigen(As); As2 <- ev$vectors %*% diag(pmax(ev$values, 1e-14)) %*% t(ev$vectors); As2 <- (As2 + t(As2)) / 2; dimnames(As2) <- list(lev, lev)
tryA("min eigenvalue 1e-14 by hand", As2)
Aext <- rbind(cbind(A0, 1), c(rep(1, 8), 0.5)); dimnames(Aext) <- list(c(lev, "E"), c(lev, "E"))
tryA("an extra level breaks PD", Aext)
