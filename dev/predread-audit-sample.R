# The audit's frmtmb.sample arm: frm_sample() draws read by frmtmb.sample and
# by brms 2.23.0, the SAME draws. Each of `nd` draws of frm_sample() is one
# chain of a Fixed_param brmsfit, so brms's post-processing runs on frmtmb's
# values. posterior_linpred, posterior_epred and fitted are compared draw by
# draw; posterior_predict by each row's mean and variance over the draws.
#   Rscript dev/predread-audit-sample.R <lib> <tag>
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r16"
tag <- if (length(args) > 1L) args[2] else "before"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r16",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
quiet <- function(e) suppressWarnings(suppressMessages(e))
ev <- function(txt, pkg) {
  eval(parse(text = txt)[[1L]], envir = list2env(list(),
       parent = if (pkg == "brms") asNamespace("brms") else
         asNamespace("frmtmb")))
}
base_data <- function(seed) {
  d <- E$eq_dat("gaussian", seed = seed)
  d$off <- round(runif(nrow(d), 0, 0.5), 3)
  d$m <- sample(0:3, nrow(d), TRUE)
  d$sx <- 0.3
  d$xo <- round(d$x + rnorm(nrow(d), 0, 0.3), 3)
  d
}
ord <- function(d) {
  d$ord <- cut(d$eta + rlogis(nrow(d)), c(-Inf, -0.5, 0.5, 1.5, Inf),
               labels = FALSE)
  d
}
shapes <- list(
  group = list(f = "bf(y ~ x + (1 + x | g))"),
  offset = list(f = "bf(y ~ x + offset(off))"),
  mo = list(f = "bf(y ~ mo(m) + x)"),
  me = list(f = "bf(y ~ me(xo, sx) + z)", noisy = TRUE),
  s = list(f = "bf(y ~ s(x, k = 5) + z)"),
  hsgp = list(f = "bf(y ~ gp(x, k = 8, c = 5/4))"),
  sigma_dpar = list(f = "bf(y ~ x, sigma ~ z)"),
  nl = list(f = "bf(y ~ a + b * x, a ~ 1 + z, b ~ 1, nl = TRUE)"),
  categorical = list(f = "bf(cat ~ x)", fam = "categorical()",
                     data = function(d) {
    d$cat <- cut(d$eta + rlogis(nrow(d)), c(-Inf, 0, 1, Inf), labels = FALSE)
    d
  }),
  cs = list(f = "bf(ord ~ x + cs(z))", fam = "sratio()", data = ord))
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
  r <- compare(a$v, b$v)
  add(shape, reader, r$ok, r$metric, r$note %||% "")
}
exact <- function(a, b) {
  d <- E$eq_rel(a, b)
  list(ok = is.finite(d) && d < 1e-8, metric = d,
       note = sprintf("dim %s vs %s", paste(dim(a), collapse = "x"),
                      paste(dim(b), collapse = "x")))
}
noisy_cmp <- function(a, b) {
  z <- (colMeans(a) - colMeans(b)) /
    sqrt(apply(a, 2, var) / nrow(a) + apply(b, 2, var) / nrow(b))
  rt <- apply(a, 2, sd) / apply(b, 2, sd)
  crit <- stats::qnorm(1 - 0.001 / (2 * ncol(a)))
  list(ok = max(abs(z)) < crit && max(abs(rt - 1)) < 5 / sqrt(2 * nrow(a)),
       metric = max(abs(rt - 1)),
       note = sprintf("max |z| %.2f, sd ratios %s", max(abs(z)),
                      paste(sprintf("%.3f", rt), collapse = " ")))
}
for (nm in names(shapes)) {
  sp <- shapes[[nm]]
  cat("shape", nm, "\n")
  d <- base_data(201)
  if (!is.null(sp$data)) d <- sp$data(d)
  nd <- data.frame(x = c(-1, 0.3, 1), z = 0.2, off = 0.25, m = 0:2, xo = 0.5,
                   sx = c(0.3, 0.6, 0.9),
                   g = factor(levels(d$g)[1], levels = levels(d$g)),
                   f = factor("a", levels = levels(d$f)))
  fam_f <- if (is.null(sp$fam)) stats::gaussian() else ev(sp$fam, "frmtmb")
  fam_b <- if (is.null(sp$fam)) stats::gaussian() else ev(sp$fam, "brms")
  form_f <- ev(sp$f, "frmtmb")
  form_b <- ev(sp$f, "brms")
  ds <- tryCatch(quiet(frm_sample(form_f, family = fam_f, data = d, chains = 1,
                                  iter = 400, seed = 7, refresh = 0)),
                 error = function(e) e)
  if (inherits(ds, "error")) {
    add(nm, "frm_sample", NA, NA_real_, conditionMessage(ds))
    next
  }
  nkeep <- if (isTRUE(sp$noisy)) 100L else 20L
  draws <- round(seq(1, posterior::ndraws(ds), length.out = nkeep))
  bfit <- tryCatch({
    prior <- E$brms_flat_prior(form_b, data = d, family = fam_b)
    code <- brms::make_stancode(form_b, data = d, family = fam_b,
                                prior = prior)
    sdat <- E$brms_standata(form_b, data = d, family = fam_b, prior = prior)
    scaffold <- quiet(brms::brm(form_b, data = d, family = fam_b, prior = prior,
                                empty = TRUE))
    idx <- frmtmb.sample:::draws_index(ds)
    inits <- lapply(draws, function(k) {
      E$stan_pars_from_fit(frmtmb.sample:::draws_fit_at(ds, k, idx), sdat, code,
                           scaffold$ranef)
    })
    mod <- E$brms_stan_model(code)
    sf <- quiet(rstan::sampling(mod, data = sdat, chains = length(inits),
                                iter = if (isTRUE(sp$noisy)) 40 else 1,
                                warmup = 0, algorithm = "Fixed_param",
                                init = inits, refresh = 0, seed = 1, cores = 1))
    scaffold$fit <- sf
    quiet(brms::rename_pars(scaffold))
  }, error = function(e) e)
  if (inherits(bfit, "error")) {
    add(nm, "brms_fixed_fit", NA, NA_real_, conditionMessage(bfit))
    next
  }
  noisy <- isTRUE(sp$noisy)
  set.seed(9)
  rp <- function(f) {
    if (!noisy) return(f())
    do.call(rbind, lapply(seq_len(40), function(i) f()))
  }
  c1 <- if (noisy) noisy_cmp else exact
  cmp(nm, "posterior_linpred(newdata)",
      function() rp(function() posterior_linpred(ds, newdata = nd,
                                                 draw_ids = draws)),
      function() brms::posterior_linpred(bfit, newdata = nd), c1)
  cmp(nm, "posterior_epred(newdata)",
      function() rp(function() posterior_epred(ds, newdata = nd,
                                               draw_ids = draws)),
      function() brms::posterior_epred(bfit, newdata = nd), c1)
  if (!noisy) {
    cmp(nm, "posterior_linpred(in sample)",
        function() posterior_linpred(ds, draw_ids = draws),
        function() brms::posterior_linpred(bfit), exact)
    cmp(nm, "posterior_epred(in sample)",
        function() posterior_epred(ds, draw_ids = draws),
        function() brms::posterior_epred(bfit), exact)
    cmp(nm, "fitted(scale = linear, summary = FALSE)",
        function() fitted(ds, newdata = nd, draw_ids = draws, scale = "linear",
                          summary = FALSE),
        function() fitted(bfit, newdata = nd, scale = "linear",
                          summary = FALSE), exact)
  }
  cmp(nm, "posterior_predict(newdata)",
      function() {
        do.call(rbind, lapply(seq_len(if (noisy) 40L else 200L), function(i) {
          posterior_predict(ds, newdata = nd, draw_ids = draws)
        }))
      },
      function() {
        do.call(rbind, lapply(seq_len(if (noisy) 1L else 200L), function(i) {
          brms::posterior_predict(bfit, newdata = nd)
        }))
      }, noisy_cmp)
  if (nm == "me") next
  cmp(nm, "conditional_effects(x, epred)",
      function() {
        ce <- conditional_effects(ds, "x", resolution = 4, ndraws = NULL,
                                  categorical = if (nm == "categorical") TRUE)
        ce[[1L]]$estimate__
      },
      function() NA, function(a, b) list(ok = NA, metric = NA_real_,
                                         note = "layout only: see test files"))
}
res <- do.call(rbind, out)
write.table(res, sprintf("dev/predread-log/audit-sample-%s.tsv", tag),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("cells:", nrow(res), "; agree:", sum(res$ok %in% TRUE), "; differ:",
    sum(res$ok %in% FALSE), "; no reading:", sum(is.na(res$ok)), "\n")
print(res[!res$ok %in% TRUE, c("shape", "reader", "metric", "note")],
      right = FALSE)
