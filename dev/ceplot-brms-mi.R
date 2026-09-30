# Lane ceplot: what brms 2.23.0 does with conditional_effects() on a
# model with mi(x, idx = ), and with one invalid effect among valid
# ones, and on an ordinal fit. brms is run with algorithm =
# "fixed_param" at its default inits, so nothing is sampled; the
# question is which calls answer and what they warn.
#   Rscript dev/ceplot-brms-mi.R > dev/ceplot-log/brms-mi.txt 2>&1
# Data seed 26 (the mi_idx_data() of test-subset-rate.R) and 2.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
options(mc.cores = 1)
try_ <- function(label, expr) {
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  if (inherits(r, "error")) {
    cat(label, ": ERROR ", conditionMessage(r), "\n", sep = "")
  } else {
    cat(label, ": OK ", paste(names(r), collapse = ", "), "\n", sep = "")
  }
  for (x in w) cat("   warning: ", x, "\n", sep = "")
  invisible(r)
}
bfix <- function(formula, data, ...) {
  suppressMessages(suppressWarnings(
    brm(formula, data = data, algorithm = "fixed_param", chains = 1,
        iter = 20, warmup = 0, init = 0, refresh = 0, seed = 1,
        silent = 2, ...)))
}

set.seed(26)
n <- 120
dm <- data.frame(g1 = sample(seq(1, n - 1, 2), n, TRUE), g2 = seq_len(n),
                 s = rep(c(TRUE, FALSE), n / 2), w = rnorm(n))
dm$x <- rnorm(n)
dm$y <- 1 + 0.5 * dm$x[match(dm$g1, dm$g2)] + 0.3 * dm$w +
  rnorm(n, sd = 0.5)
bm <- bfix(bf(y ~ mi(x, idx = g1) + w) +
             bf(x | mi() + index(g2) + subset(s) ~ 1) + set_rescor(FALSE),
           dm)
r <- try_("mi idx, default", conditional_effects(bm))
if (!inherits(r, "error")) print(head(r[[1]]))
r <- try_("mi idx, w, resp y", conditional_effects(bm, "w", resp = "y"))
if (!inherits(r, "error")) print(head(r[[1]]))
r <- try_("mi idx, x, resp y", conditional_effects(bm, "x", resp = "y"))
if (!inherits(r, "error")) print(head(r[[1]]))
r <- try_("mi idx, w, resp y, conditions g1 = g2 = 1",
          conditional_effects(bm, "w", resp = "y",
                              conditions = data.frame(g1 = 1, g2 = 1)))
if (!inherits(r, "error")) print(head(r[[1]]))

set.seed(1)
d <- data.frame(x = rnorm(100), f = factor(rep(c("a", "b"), 50)),
                z = rnorm(100))
d$y <- rnorm(100, 1 + 0.5 * d$x + (d$f == "b"))
bl <- bfix(y ~ x * f, d)
try_("effects = c(xx, x)", conditional_effects(bl, effects = c("xx", "x")))
try_("effects = c(z, x) (z in data only)",
     conditional_effects(bl, effects = c("z", "x")))
try_("effects = xx", conditional_effects(bl, effects = "xx"))
try_("effects = c(x:xx, f)", conditional_effects(bl, effects = c("x:xx", "f")))

set.seed(2)
do <- data.frame(x = rnorm(200), z = rnorm(200))
do$y <- factor(cut(do$x + rlogis(200), c(-Inf, -1, 0, 1, Inf)),
               ordered = TRUE)
bo <- bfix(y ~ x + z, do, family = cumulative())
try_("ordinal default", conditional_effects(bo))
try_("ordinal categorical = TRUE", conditional_effects(bo, categorical = TRUE))
try_("ordinal dpar = mu", conditional_effects(bo, dpar = "mu"))
try_("ordinal method = posterior_linpred",
     conditional_effects(bo, method = "posterior_linpred"))
try_("ordinal method = posterior_predict",
     conditional_effects(bo, method = "posterior_predict"))
cat("done\n")
