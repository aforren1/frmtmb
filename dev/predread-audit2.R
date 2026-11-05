# The audit's power cells, written because dev/predread-audit.R could not see 11 of
# the review's 14 mutants: no cell compared a standard error with brms, cs_cum's
# grid sat on the reference level, me had sx = 0.3 and no gr =. A record, not a
# gate: the test files gate. Same parameter values on both sides (Fixed_param at
# the frm() estimate), or the same frm_sample() draws read by both.
#   Rscript dev/predread-audit2.R <lib> <tag>
# writes dev/predread-log/audit2-<tag>.tsv. Data seed 301, sampler seed 7.
args <- commandArgs(TRUE)
lib <- args[1]
tag <- args[2]
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/wt-predread-lib",
                   "C:/Users/adf44/source/r/rellib-r16",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), "\n")
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
quiet <- function(e) suppressWarnings(suppressMessages(e))
out <- list()
add <- function(shape, cell, ok, metric, note = "") {
  out[[length(out) + 1L]] <<- data.frame(shape = shape, cell = cell, ok = ok,
    metric = signif(metric, 4), note = substr(gsub("[\r\n\t]+", " ", note), 1, 120))
}
safely <- function(shape, cell, e) {
  r <- tryCatch(quiet(e), error = function(err) {
    list(ok = FALSE, metric = NA_real_, note = paste("error:", conditionMessage(err)))
  })
  add(shape, cell, r$ok, r$metric, r$note %||% "")
}
near <- function(a, b, lim = 1e-8) {
  d <- E$eq_rel(a, b)
  list(ok = is.finite(d) && d < lim, metric = d, note = "")
}
ratio_ok <- function(a, b, lim) {
  r <- as.numeric(a) / as.numeric(b)
  list(ok = all(abs(r - 1) < lim), metric = max(abs(r - 1)),
       note = paste(sprintf("%.3f", r), collapse = " "))
}
set.seed(301)
n <- 220
d <- data.frame(x = rnorm(n), z = rnorm(n), sx = 0.5, sz = 0.3,
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(rep(sprintf("g%02d", 1:22), each = 10)))
lat <- rnorm(n)
d$xo <- lat + rnorm(n, 0, 0.5)
d$sg <- 0.3
lv <- rnorm(22)[as.integer(d$g)]
d$xg <- lv + rnorm(22, 0, 0.3)[as.integer(d$g)]
d$y <- 1 + 1.1 * lat + 0.4 * d$z + rnorm(n, 0, 0.5)
d$yg <- 1 + 1.1 * lv + 0.4 * d$z + rnorm(n, 0, 0.5)
d$ord <- cut(0.8 * d$x - 0.5 * d$z + rlogis(n), c(-Inf, -0.8, 0.2, 1.2, Inf),
             labels = FALSE)

# ---- the coefficients of a cs() term, in a standard error and on a factor
for (sp in list(list(nm = "cs", f = bf(ord ~ x + cs(z)), b = brms::bf(ord ~ x + cs(z)),
                     ff = sratio(), bb = brms::sratio(),
                     nd = data.frame(x = c(-1, 0, 1), z = c(0.5, -1, 1.5),
                                     f = factor("b", levels = levels(d$f)))),
                list(nm = "cs_cum", f = bf(ord ~ x + cs(f)),
                     b = brms::bf(ord ~ x + cs(f)), ff = acat(), bb = brms::acat(),
                     nd = data.frame(x = c(-1, 0, 1), z = 0,
                                     f = factor(c("b", "c", "b"),
                                                levels = levels(d$f)))))) {
  cat(sp$nm, "\n")
  fit <- quiet(frm(sp$f, family = sp$ff, data = d))
  bx <- quiet(E$brms_fixed_fit(sp$b, sp$bb, d, fit))
  safely(sp$nm, "frm_linpred, non-reference grid", {
    near(frm_linpred(fit, newdata = sp$nd), posterior_linpred(
      bx, newdata = sp$nd, re_formula = NA)[1, , ])
  })
  safely(sp$nm, "emmeans ~ x | rep.meas", {
    at <- list(x = c(-1, 1), z = 0.5, f = "b")
    ea <- as.data.frame(emmeans::emmeans(fit, ~ x | rep.meas, at = at))$emmean
    eb <- as.data.frame(emmeans::emmeans(bx, ~ x | rep.meas, at = at))$emmean
    near(ea, eb)
  })
  # the same draws in both packages: se.fit against brms's sd over them
  ds <- quiet(frm_sample(sp$f, family = sp$ff, data = d, chains = 1, iter = 600,
                         seed = 7, refresh = 0))
  draws <- round(seq(1, posterior::ndraws(ds), length.out = 150))
  prior <- E$brms_flat_prior(sp$b, data = d, family = sp$bb)
  code <- brms::make_stancode(sp$b, data = d, family = sp$bb, prior = prior)
  sdat <- E$brms_standata(sp$b, data = d, family = sp$bb, prior = prior)
  sc <- quiet(brms::brm(sp$b, data = d, family = sp$bb, prior = prior,
                        empty = TRUE))
  idx <- frmtmb.sample:::draws_index(ds)
  inits <- lapply(draws, function(k) {
    E$stan_pars_from_fit(frmtmb.sample:::draws_fit_at(ds, k, idx), sdat, code,
                         sc$ranef)
  })
  sc$fit <- quiet(rstan::sampling(E$brms_stan_model(code), data = sdat,
                                  chains = length(inits), iter = 1, warmup = 0,
                                  algorithm = "Fixed_param", init = inits,
                                  refresh = 0, seed = 1, cores = 1))
  bs <- quiet(brms::rename_pars(sc))
  safely(sp$nm, "se.fit against brms's sd over one set of draws", {
    a <- frm_linpred(fit, newdata = sp$nd, se.fit = TRUE)$se.fit
    b <- apply(posterior_linpred(bs, newdata = sp$nd, re_formula = NA), 2:3, sd)
    ratio_ok(a, b, 0.5)
  })
  safely(sp$nm, "frm_curve .se against brms's sd", {
    a <- frm_curve(fit, newdata = sp$nd, simultaneous = FALSE)$.se
    b <- apply(posterior_linpred(bs, newdata = sp$nd, re_formula = NA), 2:3, sd)
    ratio_ok(a, b, 0.5)
  })
}

# ---- the draw of a new row's me() value
nd <- data.frame(xo = c(-1, 0, 1), sx = c(0.4, 0.8, 1.2), z = c(0, 1, -1))
fit <- quiet(frm(bf(y ~ me(xo, sx) + z), data = d, family = gaussian()))
with_me <- function(e) {
  attach(list(me = brms::me), name = "frmtmb_a2_me", warn.conflicts = FALSE)
  on.exit(detach("frmtmb_a2_me"))
  force(e)
}
bx <- with_me(quiet(E$brms_fixed_fit(brms::bf(y ~ me(xo, sx) + z), stats::gaussian(),
                                     d, fit, ndraws = 4000)))
cat("me\n")
be <- quiet(posterior_epred(bx, newdata = nd))
safely("me", "extra_var against brms's sd at fixed values", {
  ratio_ok(sqrt(frm_lp_basis(fit, newdata = nd)$extra_var), apply(be, 2, sd), 0.06)
})
safely("me", "se.fit carries the draw and the coefficients", {
  lb <- frm_lp_basis(fit, newdata = nd)
  a <- frm_linpred(fit, newdata = nd, se.fit = TRUE)$se.fit
  ratio_ok(a, sqrt(apply(be, 2, var) + rowSums((lb$A %*% lb$V) * lb$A)), 0.06)
})
safely("me", "predict's sd against brms's", {
  set.seed(5)
  a <- predict(fit, newdata = nd, summary = FALSE, ndraws = 4000L,
               propagate_error = FALSE)
  ratio_ok(apply(a, 2, sd), apply(posterior_predict(bx, newdata = nd), 2, sd), 0.06)
})
safely("me", "rows draw apart: cor(row 1, row 2) of predict's mean", {
  set.seed(6)
  a <- t(replicate(1500, {
    sh <- fit
    sh[["me_draw"]] <- new.env(parent = emptyenv())
    as.numeric(frm_linpred(sh, newdata = nd))
  }))
  ca <- abs(cor(a[, 1], a[, 2]))
  list(ok = ca < 0.12, metric = ca, note = sprintf("brms %.3f", cor(be[, 1], be[, 2])))
})
safely("me", "emmeans mean against brms's median of its draws", {
  at <- list(z = 0, xo = 1, sx = 0.8)
  a <- as.data.frame(emmeans::emmeans(fit, ~ z, at = at))
  bb <- as.data.frame(emmeans::emmeans(bx, ~ z, at = at))
  mcse <- 1.253 * a$SE / sqrt(4000)
  zz <- abs(a$emmean - bb$emmean) / mcse
  list(ok = zz < 5, metric = zz, note = sprintf("%.4f %.4f se %.3f", a$emmean,
                                                bb$emmean, a$SE))
})
# one call read by two columns: one value, whichever column reads it
# with an interaction in the data: at a coefficient near 0 two draws read as one
d$yi <- 1 + 0.8 * lat + 0.7 * lat * d$z + rnorm(n, 0, 0.5)
fi <- quiet(frm(bf(yi ~ me(xo, sx) + me(xo, sx):z), data = d, family = gaussian()))
bi <- with_me(quiet(E$brms_fixed_fit(brms::bf(yi ~ me(xo, sx) + me(xo, sx):z),
                                     stats::gaussian(), d, fi, ndraws = 4000)))
safely("me_inter", "draw's sd against brms's", {
  ndi <- transform(nd, z = 2)
  set.seed(7)
  a <- t(replicate(2000, {
    sh <- fi
    sh[["me_draw"]] <- new.env(parent = emptyenv())
    as.numeric(frm_linpred(sh, newdata = ndi))
  }))
  ratio_ok(apply(a, 2, sd), apply(posterior_epred(bi, newdata = ndi), 2, sd), 0.08)
})
# me(gr = ): rows of one level, one value
fg <- quiet(frm(bf(yg ~ me(xg, sg, gr = g) + z), data = d, family = gaussian()))
bg <- with_me(quiet(E$brms_fixed_fit(brms::bf(yg ~ me(xg, sg, gr = g) + z),
                                     stats::gaussian(), d, fg, ndraws = 2000)))
ndg <- data.frame(xg = d$xg[c(1, 1, 11)], sg = 0.3, z = c(0, 1, 0), g = d$g[c(1, 1, 11)])
safely("me_gr", "rows of one level share a draw, with brms", {
  be <- posterior_epred(bg, newdata = ndg)
  set.seed(8)
  a <- t(replicate(1000, {
    sh <- fg
    sh[["me_draw"]] <- new.env(parent = emptyenv())
    as.numeric(frm_linpred(sh, newdata = ndg))
  }))
  ca <- cor(a[, 1], a[, 2])
  list(ok = ca > 0.999 && cor(be[, 1], be[, 2]) > 0.999, metric = 1 - ca,
       note = sprintf("%.4f %.4f", ca, cor(be[, 1], be[, 2])))
})
safely("me_gr", "extra_cov across a level", {
  E0 <- as.matrix(frm_lp_basis(fg, newdata = ndg, extra_cov = TRUE)$extra_cov)
  list(ok = isTRUE(all.equal(E0[1, 2], E0[1, 1])) && E0[1, 3] == 0, metric = E0[1, 2],
       note = "")
})

# ---- a named location's response scale, a zero standard error, and resp
fc <- quiet(frm(bf(ord ~ x + z), family = cumulative(), data = d))
bc <- quiet(E$brms_fixed_fit(brms::bf(ord ~ x + z), brms::cumulative(), d, fc))
ndc <- data.frame(x = c(-1, 0, 1), z = c(0, 0, 0))
safely("cumulative", "posterior_epred(dpar = mu) against frm_linpred(response)", {
  near(frm_linpred(fc, newdata = ndc, dpar = "mu", type = "response"),
       posterior_epred(bc, newdata = ndc, dpar = "mu")[1, ])
})
safely("cumulative", "fitted(dpar = mu) against brms", {
  near(fitted(fc, newdata = ndc, dpar = "mu")[, "Estimate"],
       posterior_epred(bc, newdata = ndc, dpar = "mu")[1, ])
})
safely("cumulative", "a curve through a row with no coefficient in it", {
  v <- frm_curve(fc, newdata = data.frame(x = c(-1, 0, 1), z = 0), dpar = "mu",
                 simultaneous = FALSE)$.estimate
  near(v, posterior_linpred(bc, newdata = ndc, dpar = "mu")[1, ])
})
safely("cumulative", "resp = character(0) as brms reads it", {
  near(frm_linpred(fc, newdata = ndc, dpar = "mu", resp = character(0)),
       posterior_linpred(bc, newdata = ndc, dpar = "mu", resp = character(0))[1, ])
})
res <- do.call(rbind, out)
write.table(res, sprintf("dev/predread-log/audit2-%s.tsv", tag), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("cells:", nrow(res), "; agree:", sum(res$ok %in% TRUE), "; differ:",
    sum(!res$ok %in% TRUE), "\n")
print(res[!res$ok %in% TRUE, ], row.names = FALSE, right = FALSE)
