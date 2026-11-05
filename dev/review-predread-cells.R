# Item 1: cells the lane's audit lacks. Same parameter values on both sides,
# offsets differing per threshold, grids with nonzero values in every term.
# Data seed 7301; tolerance 1e-8 relative on exact cells.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7301)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(sample(sprintf("g%d", 1:12), n, TRUE)))
d$m <- sample(0:3, n, TRUE)
d$lat <- d$x
d$sx <- 0.3
d$xo <- d$lat + rnorm(n, 0, 0.3)
eta <- 0.3 + 0.7 * d$x + 0.4 * d$m
d$ord <- as.integer(cut(eta + ifelse(d$z > 0, 1.2, -1.2) * rlogis(n),
                        c(-Inf, -0.5, 0.8, 2, Inf), labels = FALSE))
d$ord2 <- as.integer(cut(eta + 0.8 * as.integer(d$f) + rlogis(n),
                         c(-Inf, 0.5, 2, 3.5, Inf), labels = FALSE))
d$cat <- as.integer(cut(eta + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
d$y <- eta + rnorm(n)
d$y2 <- 0.5 * d$z + rnorm(n)
comp <- rbinom(n, 1, 0.5)
d$mord <- as.integer(cut(ifelse(comp == 1, eta + 1, -eta) + rlogis(n),
                         c(-Inf, -0.5, 0.8, 2, Inf), labels = FALSE))
nd <- data.frame(x = c(-1, 0.3, 1.2, 0.5), z = c(0.7, -1.1, 1.4, 0.2),
                 xo = c(0.5, 1, -1, 2), sx = 0.3,
                 f = factor(c("b", "c", "a", "b"), levels = levels(d$f)),
                 m = c(0L, 1L, 2L, 3L), tr = 10L,
                 g = factor("g1", levels = levels(d$g)))
cells <- list(
  cs_disc = list(f = "bf(ord ~ x + cs(z), disc ~ z)", fam = "sratio()",
                 dp = c("mu", "disc")),
  cs_two = list(f = "bf(ord2 ~ x + cs(z) + cs(f))", fam = "acat()", dp = "mu"),
  cs_mo = list(f = "bf(ord ~ mo(m) + cs(z))", fam = "acat()", dp = "mu"),
  cs_mo_in = list(f = "bf(ord ~ x + cs(mo(m)))", fam = "acat()", dp = "mu"),
  cs_re = list(f = "bf(ord ~ x + cs(z) + (1 | g))", fam = "cratio('probit')",
               dp = "mu"),
  cs_smooth = list(f = "bf(ord ~ s(x, k = 5) + cs(z))", fam = "sratio()",
                   dp = "mu"),
  cs_me = list(f = "bf(ord ~ me(xo, sx) + cs(z))", fam = "acat()", dp = "mu",
               noisy = TRUE),
  cs_cum = list(f = "bf(ord ~ x + cs(z))", fam = "cumulative()", dp = "mu"),
  hcum_cs = list(f = "bf(ord ~ x + cs(z), hu ~ z)",
                 fam = "hurdle_cumulative('probit')", dp = c("mu", "hu"),
                 data = function(d) {
                   d$ord[seq(1, n, by = 5)] <- 0L
                   d
                 }),
  mix_cs = list(f = "bf(mord ~ 1, mu1 ~ x + cs(z), mu2 ~ x)",
                fam = "mixture(sratio(), sratio())", dp = c("mu1", "mu2")),
  mix_cs2 = list(f = "bf(mord ~ 1, mu1 ~ x + cs(z), mu2 ~ x + cs(z))",
                 fam = "mixture(acat(), acat('probit'))", dp = c("mu1", "mu2")),
  cat_dp = list(f = "bf(cat ~ x, mu2 ~ x + z, mu3 ~ z)", fam = "categorical()",
                dp = c("mu2", "mu3")),
  mv_cs = list(f = "mvbf(bf(y ~ x), bf(ord ~ x + cs(z)) + sratio(), rescor = FALSE)",
               fam = "gaussian()", mv = TRUE, resp = c("y", "ord")),
  nl_cum = list(f = "bf(ord ~ a + b * x, a ~ z, b ~ 1, nl = TRUE)",
                fam = "cumulative()", dp = "mu"))
out <- list()
add <- function(cell, reader, rel, note = "") {
  out[[length(out) + 1L]] <<- data.frame(
    cell = cell, reader = reader, rel = signif(rel, 3),
    note = substr(gsub("[\r\n\t]+", " ", note), 1, 170))
}
cmp <- function(cell, reader, fa, fb) {
  a <- try_num(fa())
  b <- try_num(fb())
  note <- paste0(attr(a, "err") %||% "", " | ", attr(b, "err") %||% "")
  add(cell, reader, rd(a, b), note)
}
lay <- function(b) if (length(dim(b)) == 3L) as.numeric(b[1, , ]) else b[1, ]
for (nm in names(cells)) {
  sp <- cells[[nm]]
  cat("== ", nm, "\n")
  dd <- if (is.null(sp$data)) d else sp$data(d)
  pr <- tryCatch(rev_pair(sp$f, sp$fam, dd), error = function(e) e)
  if (inherits(pr, "error")) {
    add(nm, "pair", NA, conditionMessage(pr))
    next
  }
  fit <- pr$fit
  fx <- pr$fx
  ds <- rev_draws_at(fit, 2L)
  cat("  beta ", round(fit$estimates$beta, 3), "\n")
  if (isTRUE(sp$mv)) {
    for (r in sp$resp) {
      cmp(nm, paste0("frm_linpred ", r),
          function() frm_linpred(fit, newdata = nd, resp = r),
          function() lay(posterior_linpred(fx, newdata = nd, resp = r)))
      cmp(nm, paste0("frm_curve ", r),
          function() frm_curve(fit, newdata = nd, resp = r,
                               simultaneous = FALSE)$.estimate,
          function() lay(posterior_linpred(fx, newdata = nd, resp = r,
                                           re_formula = NA)))
      cmp(nm, paste0("frm_lp_basis ", r),
          function() frm_lp_basis(fit, newdata = nd, resp = r)$eta,
          function() lay(posterior_linpred(fx, newdata = nd, resp = r)))
      cmp(nm, paste0("fitted linear ", r), function() {
        f <- fitted(fit, newdata = nd, resp = r, scale = "linear")
        if (length(dim(f)) == 3L) f[, "Estimate", ] else f[, "Estimate"]
      }, function() lay(fitted(fx, newdata = nd, resp = r, scale = "linear",
                               summary = FALSE)))
      cmp(nm, paste0("draws linpred ", r),
          function() lay(frmtmb.sample::posterior_linpred(ds, newdata = nd,
                                                          resp = r)),
          function() lay(posterior_linpred(fx, newdata = nd, resp = r)))
      cmp(nm, paste0("draws epred ", r),
          function() lay(frmtmb.sample::posterior_epred(ds, newdata = nd,
                                                        resp = r)),
          function() lay(posterior_epred(fx, newdata = nd, resp = r)))
      cmp(nm, paste0("emmeans ~ x ", r), function() {
        as.data.frame(emmeans::emmeans(fit, ~ x, at = list(x = c(-1, 1)),
                                       resp = r))$emmean
      }, function() {
        as.data.frame(emmeans::emmeans(fx, ~ x, at = list(x = c(-1, 1)),
                                       resp = r))$emmean
      })
    }
    cmp(nm, "draws linpred mv array",
        function() lay(frmtmb.sample::posterior_linpred(ds, newdata = nd)),
        function() lay(posterior_linpred(fx, newdata = nd)))
    next
  }
  linp <- function(...) lay(posterior_linpred(fx, newdata = nd, ...))
  cmp(nm, "frm_linpred", function() frm_linpred(fit, newdata = nd, type = "link"),
      function() linp())
  cmp(nm, "frm_linpred in sample", function() frm_linpred(fit, type = "link"),
      function() lay(posterior_linpred(fx)))
  cmp(nm, "frm_lp_basis", function() frm_lp_basis(fit, newdata = nd)$eta,
      function() linp())
  cmp(nm, "frm_curve",
      function() frm_curve(fit, newdata = nd, simultaneous = FALSE)$.estimate,
      function() lay(posterior_linpred(fx, newdata = nd, re_formula = NA)))
  cmp(nm, "fitted linear", function() {
    f <- fitted(fit, newdata = nd, scale = "linear")
    if (length(dim(f)) == 3L) f[, "Estimate", ] else f[, "Estimate"]
  }, function() lay(fitted(fx, newdata = nd, scale = "linear", summary = FALSE)))
  cmp(nm, "fitted response", function() {
    f <- fitted(fit, newdata = nd)
    if (length(dim(f)) == 3L) f[, "Estimate", ] else f[, "Estimate"]
  }, function() lay(posterior_epred(fx, newdata = nd)))
  cmp(nm, "get_predict link", function() {
    marginaleffects::get_predict(fit, newdata = nd, type = "link",
                                 re_formula = NA)$estimate
  }, function() lay(posterior_linpred(fx, newdata = nd, re_formula = NA)))
  cmp(nm, "draws linpred",
      function() lay(frmtmb.sample::posterior_linpred(ds, newdata = nd)),
      function() linp())
  cmp(nm, "draws fitted linear", function() {
    lay(frmtmb.sample:::fitted.frmtmb_draws(ds, newdata = nd, scale = "linear",
                                            summary = FALSE))
  }, function() lay(fitted(fx, newdata = nd, scale = "linear", summary = FALSE)))
  cmp(nm, "draws epred",
      function() lay(frmtmb.sample::posterior_epred(ds, newdata = nd)),
      function() lay(posterior_epred(fx, newdata = nd)))
  cmp(nm, "emmeans ~ x", function() {
    as.data.frame(emmeans::emmeans(fit, ~ x, at = list(x = c(-1, 1))))$emmean
  }, function() {
    as.data.frame(emmeans::emmeans(fx, ~ x, at = list(x = c(-1, 1))))$emmean
  })
  for (dp in sp$dp) {
    cmp(nm, paste0("frm_linpred dpar ", dp),
        function() frm_linpred(fit, newdata = nd, type = "link", dpar = dp),
        function() linp(dpar = dp))
    cmp(nm, paste0("frm_linpred response dpar ", dp),
        function() frm_linpred(fit, newdata = nd, type = "response", dpar = dp),
        function() lay(posterior_epred(fx, newdata = nd, dpar = dp)))
    cmp(nm, paste0("fitted response dpar ", dp),
        function() fitted(fit, newdata = nd, dpar = dp)[, "Estimate"],
        function() lay(posterior_epred(fx, newdata = nd, dpar = dp)))
    cmp(nm, paste0("draws linpred dpar ", dp), function() {
      lay(frmtmb.sample::posterior_linpred(ds, newdata = nd, dpar = dp))
    }, function() linp(dpar = dp))
    cmp(nm, paste0("draws epred dpar ", dp), function() {
      lay(frmtmb.sample::posterior_epred(ds, newdata = nd, dpar = dp))
    }, function() lay(posterior_epred(fx, newdata = nd, dpar = dp)))
    cmp(nm, paste0("CE dpar ", dp), function() {
      conditional_effects(fit, "x", dpar = dp, resolution = 3)[[1L]]$estimate__
    }, function() {
      brms::conditional_effects(fx, "x", dpar = dp, resolution = 3)[[1L]]$estimate__
    })
    cmp(nm, paste0("emmeans dpar ", dp), function() {
      as.data.frame(emmeans::emmeans(fit, ~ x, at = list(x = c(-1, 1)),
                                     dpar = dp))$emmean
    }, function() {
      as.data.frame(emmeans::emmeans(fx, ~ x, at = list(x = c(-1, 1)),
                                     dpar = dp))$emmean
    })
  }
}
res <- do.call(rbind, out)
write.table(res, "dev/review-predread-out/cells.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("cells", nrow(res), "agree", sum(res$rel < 1e-8, na.rm = TRUE), "differ",
    sum(res$rel >= 1e-8, na.rm = TRUE), "no reading", sum(is.na(res$rel)), "\n")
print(res[is.na(res$rel) | res$rel >= 1e-8, ], right = FALSE)
