# The audit of every reader of the fitted predictor against every term kind:
# frmtmb against brms 2.23.0 at the same parameter values (a Fixed_param
# brmsfit through brms_fixed_fit(), the frm() estimate as every draw).
#   Rscript dev/predread-audit.R <lib> <tag>
# writes dev/predread-log/audit-<tag>.tsv, one row per (term kind, reader).
# ok is TRUE, FALSE, or NA where a side refuses or the reader has no analog.
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r16"
tag <- if (length(args) > 1L) args[2] else "before"
only <- if (length(args) > 2L) strsplit(args[3], ",")[[1L]]
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r16",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
quiet <- function(e) suppressWarnings(suppressMessages(e))
rel <- function(a, b) E$eq_rel(a, b)

base_data <- function(seed) {
  d <- E$eq_dat("gaussian", seed = seed)
  d$off <- round(runif(nrow(d), 0, 0.5), 3)
  d$m <- sample(0:3, nrow(d), TRUE)
  d$sx <- 0.3
  d$xo <- round(d$x + rnorm(nrow(d), 0, 0.3), 3)
  d$ng <- as.integer(d$g)
  d
}
grid_of <- function(d, n = 5) {
  data.frame(x = seq(min(d$x), max(d$x), length.out = n), z = 0.2,
             f = factor("a", levels = levels(d$f)), off = 0.25,
             m = rep_len(0:3, n), xo = 0.5, sx = 0.3,
             g = factor(levels(d$g)[1], levels = levels(d$g)), t = 1L,
             o = 1L, ntr = 10L)
}
shapes <- list(
  plain = list(f = "bf(y ~ x + z)"),
  group = list(f = "bf(y ~ x + (1 + x | g))"),
  factor = list(f = "bf(y ~ x * f)"),
  offset = list(f = "bf(y ~ x + offset(off))"),
  poly = list(f = "bf(y ~ poly(x, 2) + scale(z))"),
  mo = list(f = "bf(y ~ mo(m) + x)"),
  mo_inter = list(f = "bf(y ~ mo(m):x + z)"),
  me = list(f = "bf(y ~ me(xo, sx) + z)", noisy = TRUE),
  s = list(f = "bf(y ~ s(x, k = 5) + z)"),
  s_by = list(f = "bf(y ~ s(x, k = 5, by = f) + f)"),
  t2 = list(f = "bf(y ~ t2(x, z, k = 4))"),
  hsgp = list(f = "bf(y ~ gp(x, k = 8, c = 5/4))"),
  ar_cov = list(f = "bf(y ~ x + ar(time = t, gr = g, cov = TRUE))"),
  ar_cond = list(f = "bf(y ~ x + ar(time = t, gr = g, cov = FALSE))"),
  nl = list(f = "bf(y ~ a + b * x, a ~ 1 + z, b ~ 1, nl = TRUE)"),
  sigma_dpar = list(f = "bf(y ~ x, sigma ~ z)"),
  mixture = list(f = "bf(y ~ x + z)",
                 fam = "mixture(gaussian(), gaussian())", mix = TRUE),
  categorical = list(f = "bf(cat ~ x)", fam = "categorical()",
                     data = function(d) {
    d$cat <- cut(d$eta + rlogis(nrow(d)), c(-Inf, 0, 1, Inf), labels = FALSE)
    d
  }),
  cumulative = list(f = "bf(ord ~ x + z)", fam = "cumulative()",
                    data = function(d) {
    d$ord <- cut(d$eta + rlogis(nrow(d)), c(-Inf, -0.5, 0.5, 1.5, Inf),
                 labels = FALSE)
    d
  }),
  cs = list(f = "bf(ord ~ x + cs(z))", fam = "sratio()",
            data = function(d) {
    d$ord <- cut(d$eta + rlogis(nrow(d)), c(-Inf, -0.5, 0.5, 1.5, Inf),
                 labels = FALSE)
    d
  }),
  cs_cum = list(f = "bf(ord ~ x + cs(f))", fam = "acat()",
                data = function(d) {
    d$ord <- cut(d$eta + rlogis(nrow(d)), c(-Inf, -0.5, 0.5, 1.5, Inf),
                 labels = FALSE)
    d
  }),
  mv = list(f = "mvbf(bf(y ~ x), bf(z ~ f), rescor = FALSE)", mv = TRUE),
  mi = list(f = "bf(y ~ mi(xm) + z) + bf(xm | mi() ~ z) + set_rescor(FALSE)",
            data = function(d) {
    d$xm <- d$x
    d$xm[sample(nrow(d), 15)] <- NA
    d
  }, mv = TRUE, nogrid_mi = TRUE))
# frmtmb's bf() and brms's read one spelling, in their own package
ev <- function(txt, pkg) {
  eval(parse(text = txt)[[1L]], envir = list2env(list(),
       parent = if (pkg == "brms") asNamespace("brms") else
         asNamespace("frmtmb")))
}
out <- list()
add <- function(shape, reader, ok, metric, note = "") {
  out[[length(out) + 1L]] <<- data.frame(
    shape = shape, reader = reader, ok = ok, metric = signif(metric, 4),
    note = substr(gsub("[\r\n\t]+", " ", note), 1, 160))
}
cmp <- function(shape, reader, fa, fb, compare) {
  a <- tryCatch(list(v = quiet(fa())), error = function(e) {
    list(e = conditionMessage(e))
  })
  b <- tryCatch(list(v = quiet(fb())), error = function(e) {
    list(e = conditionMessage(e))
  })
  if (!is.null(a$e) || !is.null(b$e)) {
    add(shape, reader, NA, NA_real_, paste(
      if (!is.null(a$e)) paste("frmtmb:", a$e),
      if (!is.null(b$e)) paste("brms:", b$e)))
    return(invisible(NULL))
  }
  r <- tryCatch(compare(a$v, b$v), error = function(e) {
    list(ok = NA, metric = NA_real_, note = conditionMessage(e))
  })
  add(shape, reader, r$ok, r$metric, r$note %||% "")
}
exact <- function(a, b) {
  d <- rel(a, b)
  list(ok = is.finite(d) && d < 1e-8, metric = d,
       note = sprintf("len %d vs %d", length(a), length(b)))
}
draws_z <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  if (length(dim(b)) == 3L) b <- matrix(b, nrow = dim(b)[1L])
  cmp <- E$eq_sim_compare(a, b)
  if (is.null(cmp$crit)) return(list(ok = FALSE, metric = Inf, note = "shapes"))
  list(ok = isTRUE(cmp$ok) && isTRUE(cmp$const_agree),
       metric = max(cmp$maxz_mean, cmp$maxz_var) / cmp$crit,
       note = sprintf("z mean %.2f var %.2f crit %.2f", cmp$maxz_mean,
                      cmp$maxz_var, cmp$crit))
}
first <- function(x) if (length(dim(x)) == 3L) x[1, , ] else x[1, ]
for (nm in setdiff(names(shapes), character(0))) {
  if (!is.null(only) && !nm %in% only) next
  sp <- shapes[[nm]]
  cat("shape", nm, "\n")
  d <- base_data(101)
  if (!is.null(sp$data)) d <- sp$data(d)
  fam_f <- if (is.null(sp$fam)) stats::gaussian() else ev(sp$fam, "frmtmb")
  fam_b <- if (is.null(sp$fam)) stats::gaussian() else ev(sp$fam, "brms")
  form_f <- ev(sp$f, "frmtmb")
  form_b <- ev(sp$f, "brms")
  fit <- tryCatch(quiet(frm(form_f, family = fam_f, data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    add(nm, "frm", NA, NA_real_, conditionMessage(fit))
    next
  }
  fx <- tryCatch(quiet(E$brms_fixed_fit(form_b, fam_b, d, fit, ndraws = 4000)),
                 error = function(e) e)
  if (inherits(fx, "error")) {
    add(nm, "brms_fixed_fit", NA, NA_real_, conditionMessage(fx))
    next
  }
  nd <- grid_of(d)
  if (!is.null(sp$nogrid_mi)) nd$xm <- nd$x
  if (nm == "ar_cond" || nm == "ar_cov") {
    nd$y <- 0.1
    nd$t <- seq_len(nrow(nd))
  }
  noisy <- isTRUE(sp$noisy)
  # brms draws a new row's latent value, so frmtmb's value is the mean of that
  # draw, within its Monte Carlo error at the two-sided Bonferroni point
  first <- if (noisy) identity else function(x) {
    if (length(dim(x)) == 3L) x[1, , ] else x[1, ]
  }
  ex <- if (noisy) {
    function(a, b) {
      a <- as.numeric(a)
      z <- (a - colMeans(b)) / (apply(b, 2, sd) / sqrt(nrow(b)))
      crit <- stats::qnorm(1 - 0.001 / (2 * length(a)))
      list(ok = max(abs(z)) < crit, metric = max(abs(z)) / crit,
           note = sprintf("mean of the draw, max |z| %.2f, crit %.2f",
                          max(abs(z)), crit))
    }
  } else exact
  in_sample <- function(a, b) {
    if (noisy) return(list(ok = NA, metric = NA_real_, note = paste(
      "brms's default save_pars() keeps no latent value, so it redraws",
      "(its warning says so); frmtmb reads the fit's own")))
    exact(a, first(b))
  }
  # ---- the linear scale
  mvr <- isTRUE(sp$mv)
  cmp(nm, "frm_linpred(link, newdata)",
      function() {
        if (mvr) {
          sapply(names(fit$spec$responses), function(r) {
            as.numeric(frm_linpred(fit, newdata = nd, resp = r))
          })
        } else {
          frm_linpred(fit, newdata = nd, type = "link")
        }
      },
      function() first(brms::posterior_linpred(fx, newdata = nd,
                                               re_formula = NULL)), ex)
  cmp(nm, "frm_linpred(link, in sample)",
      function() {
        if (mvr) {
          sapply(names(fit$spec$responses), function(r) {
            as.numeric(frm_linpred(fit, resp = r))
          })
        } else {
          frm_linpred(fit, type = "link")
        }
      },
      function() brms::posterior_linpred(fx), in_sample)
  if (!mvr) {
    cmp(nm, "frm_curve(re_formula = NA)",
        function() frm_curve(fit, newdata = nd, simultaneous = FALSE)$.estimate,
        function() first(brms::posterior_linpred(fx, newdata = nd,
                                                 re_formula = NA)), ex)
    cmp(nm, "frm_lp_basis(eta)",
        function() frm_lp_basis(fit, newdata = nd)$eta,
        function() first(brms::posterior_linpred(fx, newdata = nd,
                                                 re_formula = NULL)), ex)
    cmp(nm, "frm_linpred(se.fit) vs frm_lp_basis",
        function() as.numeric(frm_linpred(fit, newdata = nd,
                                          se.fit = TRUE)$se.fit),
        function() {
          lb <- frm_lp_basis(fit, newdata = nd)
          if (is.null(lb$V_root)) {
            sqrt(pmax(rowSums((lb$A %*% lb$V) * lb$A), 0) + lb$extra_var)
          } else NA
        }, exact)
    cmp(nm, "fitted(scale = linear)",
        function() {
          f <- fitted(fit, newdata = nd, scale = "linear")
          if (length(dim(f)) == 3L) f[, "Estimate", ] else f[, "Estimate"]
        },
        function() {
          f <- fitted(fx, newdata = nd, scale = "linear", summary = FALSE)
          first(f)
        }, ex)
  }
  # ---- the response scale
  cmp(nm, "fitted(newdata)",
      function() {
        f <- fitted(fit, newdata = nd)
        if (length(dim(f)) == 3L) f[, "Estimate", ] else f[, "Estimate"]
      },
      function() first(brms::posterior_epred(fx, newdata = nd)), ex)
  cmp(nm, "fitted(in sample)",
      function() {
        f <- fitted(fit)
        if (length(dim(f)) == 3L) f[, "Estimate", ] else f[, "Estimate"]
      },
      function() brms::posterior_epred(fx), in_sample)
  if (!isTRUE(sp$fam %in% c("categorical()", "cumulative()", "sratio()",
                            "acat()")) && !mvr) {
    cmp(nm, "residuals(in sample)",
        function() residuals(fit, type = "response")[, "Estimate"],
        function() stats::residuals(fx, method = "posterior_epred",
                                    type = "ordinary", summary = FALSE),
        in_sample)
  }
  # ---- the simulators
  cmp(nm, "predict(newdata, no error)",
      function() {
        predict(fit, newdata = nd, summary = FALSE, ndraws = 4000,
                propagate_error = FALSE)
      },
      function() brms::posterior_predict(fx, newdata = nd), draws_z)
  cmp(nm, "predict(in sample, no error)",
      function() predict(fit, summary = FALSE, ndraws = 4000,
                         propagate_error = FALSE),
      function() brms::posterior_predict(fx),
      if (noisy) function(a, b) in_sample(a, b) else draws_z)
  if (!mvr) {
    cmp(nm, "simulate(newdata)",
        function() {
          s <- simulate(fit, nsim = 4000, seed = 3, newdata = nd)
          t(vapply(s, function(v) as.numeric(v), numeric(nrow(s))))
        },
        function() brms::posterior_predict(fx, newdata = nd), draws_z)
    cmp(nm, "simulate(in sample, re_formula = NULL)",
        function() {
          s <- simulate(fit, nsim = 4000, seed = 3, re_formula = NULL)
          t(vapply(s, function(v) as.numeric(v), numeric(nrow(s))))
        },
        function() brms::posterior_predict(fx),
        if (noisy) function(a, b) in_sample(a, b) else draws_z)
  }
  # ---- the interop hooks
  if (!mvr) {
    cmp(nm, "get_predict(link)",
        function() {
          marginaleffects::get_predict(fit, newdata = nd, type = "link",
                                       re_formula = NA)$estimate
        },
        function() first(brms::posterior_linpred(fx, newdata = nd,
                                                 re_formula = NA)), ex)
    cmp(nm, "get_predict(response)",
        function() {
          marginaleffects::get_predict(fit, newdata = nd,
                                       re_formula = NA)$estimate
        },
        function() first(brms::posterior_epred(fx, newdata = nd,
                                                re_formula = NA)), ex)
    evv <- if (noisy) "z" else "x"
    emat <- list(c(-1, 1))
    names(emat) <- evv
    if (nm %in% c("mo", "mo_inter")) emat$m <- 1
    emf <- stats::as.formula(paste("~", evv))
    # brms's emmean is a median over its draws; agreement within a tenth of
    # the draw's spread is all one can read for a noisy term
    near <- if (noisy) {
      function(a, b) {
        list(ok = NA, metric = max(abs(as.numeric(a) - as.numeric(b))),
             note = "median over noisy draws; read the difference")
      }
    } else exact
    cmp(nm, "emmeans ~ x (2 points)",
        function() {
          e <- as.data.frame(emmeans::emmeans(fit, emf, at = emat))
          e[[intersect(c("emmean", "response"), names(e))[1L]]]
        },
        function() {
          e <- as.data.frame(emmeans::emmeans(fx, emf, at = emat))
          e[[intersect(c("emmean", "response"), names(e))[1L]]]
        }, near)
    cmp(nm, "conditional_effects(x)",
        function() {
          ce <- conditional_effects(fit, evv, resolution = 5)
          ce[[1L]]$estimate__
        },
        function() {
          ce <- brms::conditional_effects(fx, evv, resolution = 5)
          ce[[1L]]$estimate__
        }, near)
  }
}
res <- do.call(rbind, out)
write.table(res, sprintf("dev/predread-log/audit-%s.tsv", tag), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("cells:", nrow(res), "; agree:", sum(res$ok %in% TRUE), "; differ:",
    sum(res$ok %in% FALSE), "; no reading:", sum(is.na(res$ok)), "\n")
print(res[res$ok %in% FALSE, c("shape", "reader", "metric", "note")])
print(res[is.na(res$ok), c("shape", "reader", "note")], right = FALSE)
