# Review of lane superassign: direct behavioral probes of the closures
# the lane's trace did not cover, and full-precision fit results.
# Usage: Rscript dev/review-superassign-probe.R <base|lane> <out.rds>
# base: rellib-r18 only; lane: wt-superassign-lib in front of it.
# Seeds are fixed in each probe. Output: an RDS of every result, compared
# by dev/review-superassign-probe-cmp.R.
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
out_file <- a[2]
R <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "lane") file.path(R, "wt-superassign-lib"),
            file.path(R, "rellib-r18"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.ode)
})
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "; frmtmb.ode",
    dirname(find.package("frmtmb.ode")), "\n")
res <- list()
catch <- function(expr) {
  w <- character(0)
  m <- character(0)
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                 class = "probe_error")),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")
    })
  list(value = v, warnings = w, messages = m)
}
fullp <- function(x) sprintf("%.17g", unlist(x))

# ---- P1: nlminb_trial_fn, nonfinite trials, two alive, nested ----------
tf <- frmtmb:::nlminb_trial_fn
tally <- new.env()
tally$n <- 0L
vals <- c(5, 3, NaN, Inf, NaN, 2, NA, 1, NaN)
k <- 0L
A <- tf(function(p) {
  k <<- k + 1L
  vals[k]
}, tally)
B <- tf(function(p) sum(p^2))
outA <- numeric(0)
for (i in seq_along(vals)) {
  outA <- c(outA, A$fn(c(i, -i)))
  if (i %% 2L) B$fn(c(i, 1))
}
B2 <- tf(function(p) if (p[1] > 2) NaN else p[1])
A2 <- tf(function(p) B2$fn(p) + 1)
outA2 <- vapply(1:5, function(i) A2$fn(i), 0)
res$p1 <- list(outA = outA, A = list(A$count(), A$last(), A$best()),
               B = list(B$count(), B$last(), B$best()), tally = tally$n,
               outA2 = outA2, A2 = list(A2$count(), A2$best()),
               B2 = list(B2$count(), B2$last(), B2$best()))

# ---- P2: fits whose optimizer meets nonfinite trials -------------------
set.seed(31)
n <- 200
x <- runif(n)
dg <- data.frame(x = x, y = rgamma(n, shape = 3, rate = 3 / (0.05 + 2 * x)))
p2 <- list()
for (st in list(NULL, c(0.01, 4), c(2, -1.9))) {
  f <- catch(frm(y ~ x, family = Gamma(link = "identity"), data = dg,
                 start = if (!is.null(st)) list(beta = st)))
  v <- f$value
  p2[[length(p2) + 1L]] <- if (inherits(v, "probe_error")) {
    list(error = unclass(v), w = f$warnings)
  } else {
    list(par = fullp(v$opt$par), obj = fullp(v$opt$objective),
         nonfinite = v$opt$nonfinite_trials, conv = v$opt$convergence,
         evals = v$opt$evaluations, iters = v$opt$iterations,
         w = f$warnings, m = f$messages)
  }
}
res$p2 <- p2

# ---- P3: importance objective across calls --------------------------
set.seed(7)
ng <- 25
g <- factor(rep(seq_len(ng), each = 4))
xi <- rnorm(ng * 4)
yi <- rpois(ng * 4, exp(0.3 + 0.4 * xi + rnorm(ng, 0, 0.6)[g]))
di <- data.frame(y = yi, x = xi, g = g)
p3 <- list()
for (rep_i in 1:2) {
  f <- catch(frm(y ~ x + (1 | g), family = poisson(), data = di,
                 importance = 400L))
  v <- f$value
  p3[[rep_i]] <- if (inherits(v, "probe_error")) list(error = unclass(v)) else
    list(par = fullp(v$opt$par), obj = fullp(v$opt$objective),
         ll = fullp(logLik(v)), w = f$warnings)
}
res$p3 <- p3

# ---- P4: frm_ode, a failing solve, a clean rerun, a changed one -------
blowup <- function(t, y, p) list(p[1] * y * y)
d4 <- data.frame(id = factor(c(1, 1, 2, 2)), time = c(1, 2, 1, 2))
p4 <- list()
for (pr in list(c(0.01, 0.01, 20, 20), c(0.01, 0.01, 0.02, 0.02),
                c(30, 30, 0.02, 0.02))) {
  f <- catch(frm_ode(blowup, init = list(1), times = d4$time,
                     parms = list(pr), group = d4$id, output = 1L,
                     method = "ode45"))
  p4[[length(p4) + 1L]] <- list(value = fullp(f$value), w = f$warnings,
                                failures = frm_ode_failures())
}
res$p4 <- p4

# ---- P5: frm_bind_generic, two bindings, owner unloaded and reloaded ---
fb <- frmtmb:::frm_bind_generic
e1 <- new.env()
e2 <- new.env()
fallback <- function(object, ...) "fallback"
fb(e1, "ndraws", "posterior", fallback)
fb(e2, "ndraws", "posterior", fallback)
p5 <- list()
p5$before <- c(identical(e1$ndraws, fallback), identical(e2$ndraws, fallback))
loadNamespace("posterior")
g1 <- e1$ndraws
p5$loaded <- c(identical(g1, getExportedValue("posterior", "ndraws")),
               identical(e2$ndraws, getExportedValue("posterior", "ndraws")))
p5$unload_ok <- tryCatch({unloadNamespace("posterior"); TRUE}, error = function(e) conditionMessage(e))
p5$unloaded <- c(identical(e1$ndraws, fallback), identical(e2$ndraws, fallback))
loadNamespace("posterior")
p5$reloaded <- c(identical(e1$ndraws, getExportedValue("posterior", "ndraws")),
                 identical(e1$ndraws, g1))
p5$twice <- tryCatch({
  fb(e1, "ndraws", "posterior", fallback)
  "no error"
}, error = function(e) conditionMessage(e))
res$p5 <- p5

# ---- P6: compat_rule_builder, two builders interleaved ---------------
b1 <- frmtmb::compat_rule_builder()
b2 <- frmtmb::compat_rule_builder()
b1$r("a()", "b()", "works", "one")
b2$r("c()", "d()", "refused", "two")
b1$r("e()", "f()", "untested", "three")
res$p6 <- list(b1 = b1$rules(), b2 = b2$rules(), b1again = b1$rules())

# ---- P7: me_node_value, two memos alive, calls counted -----------------
cnt <- new.env()
cnt$n <- 0L
what <- function(sh) {
  cnt$n <- cnt$n + 1L
  sh$eta_shift$by_key$a * 2
}
node1 <- frmtmb:::me_node_value(list(), list(by_key = list(a = 1:3),
                                             key = "a"), what)
node2 <- frmtmb:::me_node_value(list(), list(by_key = list(a = 10:12),
                                             key = "a"), what)
seq7 <- list(node1(0.5), node1(0.5), node2(0.5), node1(0.5, 1),
             node2(0.5), node1(c(0.5, 0.5, 0.5)), node1(0.5))
res$p7 <- list(values = seq7, calls = cnt$n)

# ---- P8: resolve_priorlist, many rows, lb/ub/assigned interplay -------
set.seed(11)
n8 <- 240
d8 <- data.frame(x = rnorm(n8), z = rnorm(n8),
                 g = factor(sample(1:12, n8, TRUE)),
                 h = factor(sample(1:6, n8, TRUE)))
d8$y <- 1 + 0.5 * d8$x - 0.3 * d8$z + rnorm(12)[d8$g] +
  rnorm(n8, 0, exp(0.2 * d8$x))
form8 <- bf(y ~ x + z + (1 + x | g) + (1 | h), sigma ~ x)
design <- frmtmb:::prior_design(form8, d8, gaussian(), list(),
                                check_response = FALSE)
pls <- list(
  all = c(set_prior("normal(0, 1)", class = "b"),
          set_prior("", class = "b", coef = "x", lb = -1, ub = 2),
          set_prior("normal(0, 3)", class = "b", coef = "x"),
          set_prior("student_t(3, 0, 1)", class = "sd"),
          set_prior("", class = "sd", group = "g", lb = 0.01),
          set_prior("normal(0, 0.5)", class = "sd", group = "g", coef = "x"),
          set_prior("", class = "sd", group = "h", ub = 5),
          set_prior("lkj(2)", class = "cor"),
          set_prior("lkj(4)", class = "cor", group = "g"),
          set_prior("normal(0, 1)", class = "b", dpar = "sigma"),
          set_prior("", class = "b", coef = "x", dpar = "sigma", lb = -3),
          set_prior("normal(1, 2)", class = "Intercept"),
          set_prior("normal(0, 2)", class = "Intercept", dpar = "sigma")),
  order1 = c(set_prior("", class = "sd", lb = 0.1, ub = 3),
             set_prior("exponential(1)", class = "sd"),
             set_prior("", class = "sd", group = "g", coef = "Intercept",
                       ub = 2)),
  order2 = c(set_prior("exponential(1)", class = "sd"),
             set_prior("", class = "sd", lb = 0.1, ub = 3),
             set_prior("gamma(2, 1)", class = "sd", group = "g",
                       coef = "Intercept"),
             set_prior("", class = "sd", group = "g", coef = "Intercept",
                       ub = 2)),
  bounds_only = c(set_prior("", class = "b", lb = -5, ub = 5),
                  set_prior("", class = "sd", ub = 10),
                  set_prior("", class = "b", coef = "z", ub = 0)),
  bad_group = c(set_prior("normal(0, 1)", class = "sd", group = "nope")),
  sd_lb0 = c(set_prior("normal(0, 1)", class = "sd", lb = 0),
             set_prior("", class = "sd", group = "h", lb = 0.5)))
p8 <- list()
for (nm in names(pls)) {
  r <- catch(frmtmb:::resolve_prior_input(design, pls[[nm]]))
  v <- r$value
  p8[[nm]] <- if (inherits(v, "probe_error")) list(error = unclass(v)) else
    list(entries = v$entries, lower = fullp(v$lower), lower_n = names(v$lower),
         upper = fullp(v$upper), upper_n = names(v$upper), w = r$warnings)
}
f8 <- catch(frm(form8, data = d8, prior = pls$all))
p8$fit <- if (inherits(f8$value, "probe_error")) unclass(f8$value) else
  list(par = fullp(f8$value$opt$par), obj = fullp(f8$value$opt$objective),
       w = f8$warnings)
res$p8 <- p8

# ---- P9: full precision on optimizer-heavy fits ------------------------
set.seed(5)
n9 <- 300
d9 <- data.frame(x = runif(n9, 0, 3), g = factor(sample(1:15, n9, TRUE)),
                 m = factor(sample(c("lo", "mid", "hi"), n9, TRUE),
                            levels = c("lo", "mid", "hi"), ordered = TRUE))
d9$y <- 2 * exp(0.5 * d9$x) + rnorm(15, 0, 0.5)[d9$g] + rnorm(n9, 0, 0.4)
d9$yc <- rpois(n9, exp(0.2 + 0.3 * as.integer(d9$m)))
fits <- list(
  mixed = quote(frm(y ~ x + (1 + x | g), data = d9)),
  nl = quote(frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
                 data = d9, start = list(beta = c(1, 0.1)))),
  mo = quote(frm(yc ~ mo(m), family = poisson(), data = d9)),
  sigma = quote(frm(bf(y ~ x, sigma ~ x), data = d9)))
p9 <- list()
for (nm in names(fits)) {
  f <- catch(eval(fits[[nm]]))
  v <- f$value
  p9[[nm]] <- if (inherits(v, "probe_error")) list(error = unclass(v)) else
    list(par = fullp(v$opt$par), obj = fullp(v$opt$objective),
         evals = v$opt$evaluations, mo = v$opt$mo_search,
         nonfinite = v$opt$nonfinite_trials, w = f$warnings,
         vc = tryCatch(utils::capture.output(print(VarCorr(v), digits = 17)),
                       error = function(e) conditionMessage(e)))
}
res$p9 <- p9

saveRDS(res, out_file)
cat("RESULT probes written:", paste(names(res), collapse = ", "), "\n")
