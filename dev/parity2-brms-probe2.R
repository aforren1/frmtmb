# Lane parity2: the sections of dev/parity2-brms-probe.R after 3919, rerun
# Lane parity2: what brms 2.23.0 does on each item of the lane, run on the
# data of dev/parity2-data.R. Fits that only need names or layouts use
# algorithm = "fixed_param" at init 0 (nothing is sampled); the new-level
# band (1497) and the loo subset (5331) use a sampled fit, 2 chains of
# 1000 iterations, seed 1.
#   Rscript dev/parity2-brms-probe2.R > dev/parity2-out/brms-probe2.txt 2>&1
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
options(mc.cores = 1)
cat("brms", as.character(packageVersion("brms")), "\n")
here <- "C:/Users/adf44/source/r/frmtmb-wt-parity2/dev"
source(file.path(here, "parity2-data.R"))
try_ <- function(label, expr) {
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  if (inherits(r, "error")) cat(label, ": ERROR ", conditionMessage(r), "\n",
                                sep = "")
  for (x in w) cat("   warning: ", substr(x, 1, 200), "\n", sep = "")
  invisible(r)
}
bfix <- function(formula, data, ...) {
  suppressMessages(suppressWarnings(
    brm(formula, data = data, algorithm = "fixed_param", chains = 1,
        iter = 20, warmup = 0, init = 0, refresh = 0, seed = 1,
        silent = 2, ...)))
}
ce_names <- function(label, fit, ...) {
  r <- try_(label, conditional_effects(fit, ...))
  if (!inherits(r, "error")) {
    cat(label, ": ", paste(names(r), collapse = " "), "\n", sep = "")
  }
  invisible(r)
}

cat("\n#### 4160 attr(x, 'min') of a mo() variable\n")
set.seed(9)
d3 <- data.frame(k = rep(1:4, 10))
d3$y <- 0.5 * d3$k + rnorm(40)
attr(d3$k, "min") <- 0
sd3 <- make_standata(y ~ mo(k), d3)
cat("fit-time Jmo:", sd3$Jmo, " Xmo_1 head:", head(sd3$Xmo_1, 8), "\n")
f3 <- try_("mo fit attr min 0", bfix(y ~ mo(k), d3))
if (!inherits(f3, "error")) {
  nd <- data.frame(k = 1:4, y = 0)
  cat("attr-min fit, newdata 1..4 Xmo_1:", standata(f3, newdata = nd)$Xmo_1,
      " Jmo", standata(f3, newdata = nd)$Jmo, "\n")
  nd0 <- data.frame(k = 0:4, y = 0)
  r <- try_("newdata 0..4", standata(f3, newdata = nd0))
  if (!inherits(r, "error")) cat("newdata 0..4 Xmo_1:", r$Xmo_1, "\n")
}
d4 <- data.frame(k = rep(2:5, 10), y = rnorm(40))
sd4 <- make_standata(y ~ mo(k), d4)
cat("no attr, 2..5: Jmo", sd4$Jmo, " Xmo_1 head:", head(sd4$Xmo_1, 4), "\n")
d5 <- d4
attr(d5$k, "min") <- 3
r <- try_("attr min above the data", make_standata(y ~ mo(k), d5))
if (!inherits(r, "error")) cat("attr min 3 on 2..5: Xmo_1 head",
                               head(r$Xmo_1, 4), " Jmo", r$Jmo, "\n")

cat("\n#### 1560 lf() on a dpar bf() already sets\n")
f15 <- try_("lf", bf(y ~ x, sigma ~ z) + lf(sigma ~ w))
if (!inherits(f15, "error")) {
  print(f15)
  cat("pforms sigma:", deparse(f15$pforms$sigma), "\n")
}

cat("\n#### 5337 drop_unused_levels = FALSE\n")
sdu <- make_standata(y ~ f + (1 | g), du, drop_unused_levels = FALSE)
cat("N_1", sdu$N_1, " K", sdu$K, " X cols", colnames(sdu$X), "\n")
cat("X col sums", colSums(sdu$X), "\n")
sdu2 <- make_standata(y ~ f + (1 + f | g), du, drop_unused_levels = FALSE)
cat("slope: N_1", sdu2$N_1, " M_1", sdu2$M_1, " K", sdu2$K, "\n")
sdu3 <- make_standata(y ~ f + (1 | g), du)
cat("default: N_1", sdu3$N_1, " K", sdu3$K, "\n")
fdu <- try_("dul fit", bfix(y ~ f + (1 | g), du, drop_unused_levels = FALSE))
if (!inherits(fdu, "error")) {
  cat("variables:", paste(grep("^(b_|r_g|sd_)", variables(fdu), value = TRUE),
                          collapse = " "), "\n")
  print(dimnames(ranef(fdu)$g))
}

cat("\n#### 5315 b_Intercept order under 0 + Intercept\n")
fi <- try_("icpt", bfix(bf(y ~ 0 + Intercept + x, sigma ~ z), d))
if (!inherits(fi, "error")) cat("variables:", paste(variables(fi),
                                                   collapse = " "), "\n")
fi2 <- try_("icpt mu only", bfix(y ~ 0 + Intercept + x, d))
if (!inherits(fi2, "error")) cat("variables:", paste(variables(fi2),
                                                    collapse = " "), "\n")
fi3 <- try_("icpt sigma 0+I", bfix(bf(y ~ x, sigma ~ 0 + Intercept + z), d))
if (!inherits(fi3, "error")) cat("variables:", paste(variables(fi3),
                                                    collapse = " "), "\n")

cat("\n#### 5325 ordinal posterior_epred draw margin\n")
fo <- try_("cum", bfix(yo ~ x, dm, family = cumulative()))
if (!inherits(fo, "error")) {
  pe <- posterior_epred(fo)
  cat("dim", dim(pe), "\n")
  str(dimnames(pe))
  pe2 <- posterior_epred(fo, newdata = dm[1:3, ])
  str(dimnames(pe2))
}
fc <- try_("cat", bfix(yo ~ x, dm, family = categorical()))
if (!inherits(fc, "error")) str(dimnames(posterior_epred(fc)))

cat("\n#### 5323 insight::get_predicted on multinomial\n")
fmn <- try_("multinomial", bfix(Y | trials(size) ~ x, dm,
                                family = multinomial()))
if (!inherits(fmn, "error")) {
  gp <- try_("get_predicted", insight::get_predicted(fmn))
  if (!inherits(gp, "error")) {
    cat("class", class(gp), " length", length(gp), "\n")
    str(gp)
  }
  gp2 <- try_("get_predicted epred",
              insight::get_predicted(fmn, predict = "expectation"))
  if (!inherits(gp2, "error")) str(gp2)
}

cat("\n#### 5300 latent residual err[k]\n")
dp <- data.frame(t = rep(1:8, 6), g = factor(rep(1:6, each = 8)),
                 x = dg$x)
set.seed(4)
dp$y <- rpois(48, exp(0.5 + 0.3 * dp$x))
fa <- try_("ar cov=FALSE", bfix(y ~ x + ar(time = t, gr = g, cov = FALSE),
                                dp, family = poisson()))
if (!inherits(fa, "error")) {
  v <- variables(fa)
  cat("err/zerr/sderr/ar:", paste(grep("^(err|zerr|sderr|ar)", v,
                                       value = TRUE)[1:6], collapse = " "),
      " count err", sum(grepl("^err\\[", v)), " zerr",
      sum(grepl("^zerr\\[", v)), "\n")
  dr <- as_draws_df(fa)
  cat("err[1..3] - sderr*zerr[1..3] at draw 1:",
      if (all(c("err[1]", "zerr[1]") %in% names(dr))) {
        unlist(dr[1, c("err[1]", "err[2]", "err[3]")]) -
          dr$sderr[1] * unlist(dr[1, c("zerr[1]", "zerr[2]", "zerr[3]")])
      } else "n/a", "\n")
}
fa2 <- try_("ar cov=TRUE", bfix(y ~ x + ar(time = t, gr = g, cov = TRUE),
                                dp, family = poisson()))
if (!inherits(fa2, "error")) {
  v <- variables(fa2)
  cat("cov=TRUE err count", sum(grepl("^err\\[", v)), " zerr count",
      sum(grepl("^zerr\\[", v)), "\n")
}

cat("\n#### 2980 brm_multiple + int_conditions\n")
fm <- try_("brm_multiple", suppressMessages(suppressWarnings(
  brm_multiple(y ~ x * z, data = list(d, d), algorithm = "fixed_param",
               chains = 1, iter = 20, warmup = 0, init = 0, refresh = 0,
               seed = 1, silent = 2))))
if (!inherits(fm, "error")) {
  r <- ce_names("multiple ce int_conditions", fm, effects = "x:z",
                int_conditions = list(z = c(-1, 1)))
  if (!inherits(r, "error")) print(unique(r[[1]]$z))
}

cat("\n#### 1497 new levels: sampled fit y ~ x + (1 | g)\n")
fs1 <- suppressMessages(suppressWarnings(
  brm(y ~ x + (1 | g), data = dg, chains = 2, iter = 1000, refresh = 0,
      seed = 1, silent = 2)))
for (snl in c("uncertainty", "gaussian")) {
  set.seed(2)
  ce <- conditional_effects(fs1, "x", re_formula = NULL,
                            sample_new_levels = snl,
                            int_conditions = list(x = c(-1, 0, 1)))
  cat(snl, "width at x = -1, 0, 1:",
      format(ce$x$upper__ - ce$x$lower__, digits = 4),
      " estimate:", format(ce$x$estimate__, digits = 4), "\n")
}
set.seed(2)
ce0 <- conditional_effects(fs1, "x", re_formula = NULL,
                           int_conditions = list(x = c(-1, 0, 1)))
cat("default width:", format(ce0$x$upper__ - ce0$x$lower__, digits = 4),
    "\n")
cat("sd_g mean:", mean(as_draws_df(fs1)$sd_g__Intercept), "\n")

cat("\n#### 5331 loo_*() on a subset of the fitted rows\n")
sub <- dg[1:10, ]
for (fn in c("loo_predict", "loo_epred", "loo_linpred",
             "loo_predictive_interval")) {
  a <- try_(paste(fn, "full"), suppressWarnings(get(fn)(fs1)))
  b <- try_(paste(fn, "subset"), suppressWarnings(get(fn)(fs1,
                                                          newdata = sub)))
  if (!inherits(a, "error") && !inherits(b, "error")) {
    ac <- if (is.matrix(a)) a[1:10, , drop = FALSE] else a[1:10]
    cat(fn, ": subset class", class(b)[1], " dim/len",
        if (is.matrix(b)) dim(b) else length(b),
        " max |subset - full[1:10]|", max(abs(b - ac)), "\n")
  }
}
r2a <- try_("loo_R2 full", suppressWarnings(loo_R2(fs1)))
r2b <- try_("loo_R2 subset", suppressWarnings(loo_R2(fs1, newdata = sub)))
if (!inherits(r2b, "error")) print(r2b)
cat("done\n")
