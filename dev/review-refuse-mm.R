# Item 4: mm() with |ID| and cov against brms. Unbalanced groups, a level only
# in the second member, different weights per term, A with extra or missing
# levels, a non-PD A, by, new data, prior tables. Seed 9800.
source("dev/review-refuse-lib.R")
set.seed(9800)
n <- 150
lev <- paste0("L", sprintf("%02d", 1:12))
# unbalanced: L11 and L12 only as a second member; L01 dominates the first
g1 <- sample(lev[1:10], n, TRUE, prob = c(0.35, rep(0.65 / 9, 9)))
g2 <- sample(lev[2:12], n, TRUE, prob = c(rep(0.6 / 9, 9), 0.2, 0.2))
d <- data.frame(g1 = g1, g2 = g2, x = rnorm(n), w = rnorm(n),
                w1 = runif(n, 0.2, 1), w2 = runif(n, 0.2, 1),
                v1 = runif(n, 0.2, 1), v2 = runif(n, 0.2, 1),
                x1 = rnorm(n), x2 = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)))
re <- rnorm(12, 0, 0.8); names(re) <- lev; res <- rnorm(12, 0, 0.3); names(res) <- lev
d$y <- 1 + 0.5 * d$x + (d$w1 * re[d$g1] + d$w2 * re[d$g2]) / (d$w1 + d$w2) +
  rnorm(n, 0, exp(-0.2 + (res[d$g1] + res[d$g2]) / 2))
show <- function(lab, r) {
  if (is_err(r)) { pr("%-40s ERROR %s", lab, msg(r, 170)); return(invisible()) }
  pr("%-40s verdict=%s const=%s pts=[%s] grad=%s errs: frm[%s] brms[%s] l1[%s]", lab, r$verdict,
     format(r$l1_const %||% NA, digits = 3),
     paste(format(unlist(r$l1_extra_const) %||% NA, digits = 3), collapse = ","),
     format(r$l1_grad %||% NA, digits = 3),
     substr(r$frm_error %||% "", 1, 110), substr(r$brms_error %||% "", 1, 110),
     substr(r$l1_error %||% "", 1, 110))
}
sp <- function(bf, data = d, ...) c(list(bf = bf, data = data), list(...))

# levels in frmtmb's pooled block: first member's then second's, or sorted?
A0 <- exp(-abs(outer(1:12, 1:12, "-")) / 3); dimnames(A0) <- list(lev, lev)
Ax <- rbind(cbind(A0, 0.1), 0.1); Ax[13, 13] <- 1
dimnames(Ax) <- list(c(lev, "EXTRA"), c(lev, "EXTRA"))
Ax <- Ax[c(13, 12:1), c(13, 12:1)]   # extra level, reordered

show("|ID| across dpars, unbalanced", tr(meas(sp(
  "bf(y ~ x + (1 + x | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2)))"))))
show("|ID| different weights per term", tr(meas(sp(
  "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(v1, v2))))"))))
show("|ID| scale = FALSE", tr(meas(sp(
  "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2), scale = FALSE)), sigma ~ 1 + (1 | p | mm(g1, g2)))"))))
show("id = spelling two terms in mu", tr(meas(sp(
  "bf(y ~ x + (1 | mm(g1, g2, id = 'p')) + (0 + w | mm(g1, g2, id = 'p')))"))))
show("|ID| two responses", tr(meas(c(sp(
  "mvbf(bf(y ~ x + (1 | p | mm(g1, g2))), bf(w ~ x + (1 | p | mm(g1, g2))), rescor = FALSE)"), list(family = NA)))))
show("mm cov = A (reordered, extra level)", tr(meas(sp(
  "bf(y ~ x + (1 | mm(g1, g2, cov = A)))", data2 = list(A = Ax)))))
show("mm cov = A two columns", tr(meas(sp(
  "bf(y ~ x + (1 + x | mm(g1, g2, cov = A)))", data2 = list(A = A0)))))
show("mm cov = A weights, mmc", tr(meas(sp(
  "bf(y ~ x + (1 + mmc(x1, x2) | mm(g1, g2, weights = cbind(w1, w2), cov = A)))", data2 = list(A = A0)))))
show("mm cov = A with ||", tr(meas(sp(
  "bf(y ~ x + (1 + x || mm(g1, g2, cov = A)))", data2 = list(A = A0)))))
show("mm cov = A cor = FALSE", tr(meas(sp(
  "bf(y ~ x + (1 + x | mm(g1, g2, cor = FALSE, cov = A)))", data2 = list(A = A0)))))
show("mm id + cov", tr(meas(sp(
  "bf(y ~ x + (1 | p | mm(g1, g2, cov = A)), sigma ~ 1 + (1 | p | mm(g1, g2, cov = A)))", data2 = list(A = A0)))))
show("|ID| with by", tr(meas(sp(
  "bf(y ~ x + (1 | p | mm(g1, g2, by = f)), sigma ~ 1 + (1 | p | mm(g1, g2, by = f)))"))))
show("mm by without id, control", tr(meas(sp(
  "bf(y ~ x + (1 | mm(g1, g2, by = f)))"))))

# refusals and the absent case, both packages
fr <- function(lab, form, data = d, data2 = NULL) {
  a <- tr(frm(eval(parse(text = form), envir = asNamespace("frmtmb")), data = data, data2 = data2))
  bf0 <- eval(parse(text = form), envir = asNamespace("brms"))
  b <- tr(brms::make_standata(bf0, data = data, data2 = data2))
  pr("%-40s frm %s | brms %s", lab, msg(a, 130), msg(b, 130))
  invisible(list(a, b))
}
Am <- A0[-3, -3]            # a level missing from A
Anp <- A0; Anp[1, 2] <- Anp[2, 1] <- 2   # symmetric, not PD
Aas <- A0; Aas[1, 2] <- 0.9              # not symmetric
Ase <- A0; Ase[1, ] <- Ase[2, ]; Ase[, 1] <- Ase[, 2]   # singular
Anon <- unname(A0)
fr("A missing a level", "bf(y ~ x + (1 | mm(g1, g2, cov = A)))", data2 = list(A = Am))
fr("A not PD", "bf(y ~ x + (1 | mm(g1, g2, cov = A)))", data2 = list(A = Anp))
fr("A not symmetric", "bf(y ~ x + (1 | mm(g1, g2, cov = A)))", data2 = list(A = Aas))
fr("A singular", "bf(y ~ x + (1 | mm(g1, g2, cov = A)))", data2 = list(A = Ase))
fr("A with no dimnames", "bf(y ~ x + (1 | mm(g1, g2, cov = A)))", data2 = list(A = Anon))
fr("control: gr(g1, cov = Anp)", "bf(y ~ x + (1 | gr(g1, cov = A)))", data2 = list(A = Anp))
fr("mm by + cov", "bf(y ~ x + (1 | mm(g1, g2, by = f, cov = A)))", data2 = list(A = A0))
fr("a key over other members", "bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g2, g1)))")
fr("a key on mm and on a factor", "bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ 1 + (1 | p | g1))")
fr("key twice", "bf(y ~ x + (1 | p | mm(g1, g2, id = 'p')))")
fr("id = NA", "bf(y ~ x + (1 | mm(g1, g2, id = NA)))")
fr("two different cov under one key", "bf(y ~ x + (1 | p | mm(g1, g2, cov = A)), sigma ~ 1 + (1 | p | mm(g1, g2)))", data2 = list(A = A0))
fr("mm dist = student", "bf(y ~ x + (1 | mm(g1, g2, dist = 'student')))")
fr("mm pw", "bf(y ~ x + (1 | mm(g1, g2, pw = rep(1, 12))))")

# new data: the merged block and A, by hand from the fit's own effects
fit <- frm(bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))),
              sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(v1, v2)))), data = d)
re1 <- ranef(fit)
print(names(re1)); str(dimnames(re1[[1]]))
nd <- d[1:6, ]
ft <- fitted(fit, newdata = nd, dpar = "mu")
own <- fitted(fit, dpar = "mu")[1:6, 1]
pr("id, weights: fitted(newdata) vs fitted at the row: %.3e", max(abs(ft[, 1] - own)))
fs <- tr(fitted(fit, newdata = nd, dpar = "sigma"))
own <- tr(fitted(fit, dpar = "sigma"))
pr("id, weights: sigma fitted(newdata) vs own: %s", if (is_err(fs) || is_err(own)) paste(msg(fs), msg(own)) else format(max(abs(fs[1:6, 1] - own[1:6, 1])), digits = 3))
# a hand mean: scale the weights, take the block's mu column
cf <- fixef(fit)
print(round(cf[, 1], 3))
rmat <- re1[[1]]
nmz <- dimnames(rmat)
mu_col <- grep("mu|Intercept", nmz[[3]], value = TRUE)
print(nmz[[3]])
one_col <- nmz[[3]][1]
rr <- rmat[, "Estimate", one_col]
hand <- cf["Intercept", 1] + cf["x", 1] * nd$x +
  (nd$w1 * rr[nd$g1] + nd$w2 * rr[nd$g2]) / (nd$w1 + nd$w2)
pr("hand mu vs fitted(newdata): %.3e", max(abs(hand - ft[, 1])))
nd2 <- nd; nd2$g2[1] <- "NEW"; nd2$g1[2] <- "NEW"; nd2$g1[3] <- "NEW2"; nd2$g2[3] <- "NEW3"
a <- tr(fitted(fit, newdata = nd2, dpar = "mu"))
pr("new level without allow_new_levels: %s", msg(a, 110))
a <- tr(fitted(fit, newdata = nd2, dpar = "mu", allow_new_levels = TRUE))
pr("allow_new_levels: %s; row 1 vs hand %s", msg(a, 110), if (is_err(a)) "-" else format(a[1, 1] - (cf["Intercept", 1] + cf["x", 1] * nd2$x[1] + nd2$w1[1] * rr[nd2$g1[1]] / (nd2$w1[1] + nd2$w2[1])), digits = 3))
for (snl in c("uncertainty", "gaussian", "old_levels")) {
  a <- tr(fitted(fit, newdata = nd2, dpar = "mu", allow_new_levels = TRUE, sample_new_levels = snl))
  pr("  sample_new_levels = %s: %s %s", snl, msg(a, 110), if (is_err(a)) "" else paste(round(a[1:3, 1], 3), collapse = " "))
}
a <- tr(fitted(fit, newdata = nd, re_formula = NA, dpar = "mu"))
pr("re_formula = NA vs fixed: %s %s", msg(a, 100), if (is_err(a)) "" else format(max(abs(a[, 1] - (cf["Intercept", 1] + cf["x", 1] * nd$x))), digits = 3))
a <- tr(fitted(fit, newdata = nd, re_formula = ~ (1 | p | mm(g1, g2, weights = cbind(w1, w2))), dpar = "mu"))
pr("re_formula written back: %s %s", msg(a, 150), if (is_err(a)) "" else format(max(abs(a[, 1] - fitted(fit, newdata = nd, dpar = "mu")[, 1])), digits = 3))
# a newdata weight changes the row, per term
nd3 <- nd; nd3$w1 <- 1; nd3$w2 <- 0; nd3$v1 <- 0; nd3$v2 <- 1
a <- fitted(fit, newdata = nd3, dpar = "mu")
pr("newdata weights: mu vs hand %.3e", max(abs(a[, 1] - (cf["Intercept", 1] + cf["x", 1] * nd3$x + rr[nd3$g1]))))
a <- tr(predict(fit, newdata = nd2, allow_new_levels = TRUE))
pr("predict(new levels) under id: %s", msg(a, 130))

fa <- frm(bf(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2), cov = A))), data = d, data2 = list(A = Ax))
ra <- ranef(fa)[[1]]
print(dimnames(ra)[[1]])
pr("levels of the block sorted as A? %s", identical(dimnames(ra)[[1]], rownames(Ax)))
# is the pooled block's A the right sub-block by name? density by hand
pv <- fa$obj$env$last.par.best
lv <- fa$frame$re_blocks[[1]]$levels %||% dimnames(ra)[[1]]
sig <- VarCorr(fa)[[1]]$sd[1]
b <- ra[, "Estimate", 1]
Asub <- Ax[lv, lv]
pr("A sub-block has dim %d; levels %s", nrow(Asub), paste(head(lv, 3), collapse = ","))
a <- tr(fitted(fa, newdata = d[1:4, ]))
pr("fitted(newdata), mm cov: %s, vs own %s", msg(a, 90), if (is_err(a)) "-" else format(max(abs(a[, 1] - fitted(fa)[1:4, 1])), digits = 3))
nd4 <- d[1:4, ]; nd4$g2[1] <- "EXTRA"; nd4$g1[2] <- "NEWZ"
a <- tr(fitted(fa, newdata = nd4, allow_new_levels = TRUE))
pr("cov, a new level that A names and one it does not: %s %s", msg(a, 130), if (is_err(a)) "" else paste(round(a[1:2, 1], 3), collapse = " "))

# prior tables and names
forms <- list(
  id = "bf(y ~ x + (1 + x | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2)))",
  idw = "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(v1, v2))))",
  cov = "bf(y ~ x + (1 + x | mm(g1, g2, cov = A)))",
  mmc = "bf(y ~ x + (1 + mmc(x1, x2) | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2)))",
  by = "bf(y ~ x + (1 | p | mm(g1, g2, by = f)), sigma ~ 1 + (1 | p | mm(g1, g2, by = f)))")
key <- function(t) {
  t <- as.data.frame(t)
  paste(t$class, t$coef, t$group, t$resp, t$dpar, t$nlpar, sep = "|")
}
for (nm in names(forms)) {
  a <- tr(default_prior(eval(parse(text = forms[[nm]]), envir = asNamespace("frmtmb")), data = d, data2 = list(A = A0)))
  b <- tr(brms::default_prior(eval(parse(text = forms[[nm]]), envir = asNamespace("brms")), data = d, data2 = list(A = A0)))
  if (is_err(a) || is_err(b)) { pr("table %s: frm %s brms %s", nm, msg(a), msg(b)); next }
  ka <- key(a); kb <- key(b)
  pr("table %s: same set %s same order %s; only frm [%s] only brms [%s]", nm, setequal(ka, kb), identical(ka, kb),
     paste(setdiff(ka, kb), collapse = "; "), paste(setdiff(kb, ka), collapse = "; "))
}
