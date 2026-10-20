# Reviewer: brms 2.23.0 on refcat = "b" with a group term, 3 and 4
# categories, a multinomial with varying trials, and make_stan_names
# edge cases. Log dev/catmulti-rev-out/brms.txt; fits saved as rds.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
out <- "dev/catmulti-rev-out"
source("dev/catmulti-rev-data.R")
sec <- function(t) cat("\n========", t, "========\n")
shape <- function(nm, x) {
  cat(nm, ": dim", paste(dim(x), collapse = " x "), "\n")
  dn <- dimnames(x)
  for (k in seq_along(dn)) {
    cat("  dimnames[[", k, "]]: ",
        if (is.null(dn[[k]])) "NULL" else
          paste(utils::head(dn[[k]], 8), collapse = ", "), "\n", sep = "")
  }
}
try_shape <- function(nm, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r) && length(r) == 1L) cat(nm, ":", r, "\n") else shape(nm, r)
  invisible(r)
}
ctl <- list(chains = 2, iter = 1500, warmup = 750, refresh = 0,
            backend = "rstan", seed = 7, cores = 2)

## --- names only (no compile) -----------------------------------------
sec("make_stan_names edge cases")
for (lv in list(c("a.b", "ab", "c"), c("a b", "1x", "c"), c("a_b", "c.d", "e"),
                c("x", "x.1", "x1"))) {
  d <- data.frame(y = factor(rep(lv, 10), levels = lv), x = rnorm(30))
  r <- tryCatch(unique(default_prior(y ~ x, data = d,
                                     family = categorical())$dpar),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("levels", paste(lv, collapse = " | "), "->", paste(r, collapse = " "), "\n")
  r <- tryCatch(unique(default_prior(y ~ x, data = d,
                                     family = categorical(refcat = NA))$dpar),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("  refcat NA ->", paste(r, collapse = " "), "\n")
}
dm <- rev_mult_data(31)
for (f in c("cbind(y1, y2, y3) | trials(N) ~ x",
            "cbind(y1, y2 + 0, y3) | trials(N) ~ x",
            "cbind(a = y1, `b c` = y2, y3) | trials(N) ~ x",
            "Yc | trials(N) ~ x")) {
  r <- tryCatch(unique(default_prior(as.formula(f), data = dm,
                                     family = multinomial())$dpar),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("multinomial", f, "->", paste(r, collapse = " "), "\n")
}
r <- tryCatch(multinomial(3), error = function(e) paste("ERROR:", conditionMessage(e)))
cat("multinomial(3):", if (is.character(r)) r else paste(class(r), r$link), "\n")
r <- tryCatch(categorical(refcat = 2)$refcat, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("categorical(refcat = 2)$refcat:", r, "\n")
dm2 <- dm; colnames(dm2$Yc) <- c("1", "2", "3")
r <- tryCatch(unique(default_prior(Yc | trials(N) ~ x, data = dm2,
                                   family = multinomial(refcat = 2))$dpar),
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat("multinomial(refcat = 2) on columns 1..3 ->", paste(r, collapse = " "), "\n")

## --- logistic_normal mvbf stancode, for the code_fix check ------------
sec("mvbf logistic_normal stancode")
set.seed(46)
n <- 50
dl <- data.frame(x = rnorm(n))
E <- cbind(0, 0.4 + 0.5 * dl$x + rnorm(n), -0.3 + rnorm(n))
Y <- exp(E) / rowSums(exp(E)); colnames(Y) <- c("a", "b", "c")
dl$Y <- Y; dl$w <- rnorm(n)
sc <- stancode(bf(Y ~ x, family = logistic_normal()) +
                 bf(w ~ x, family = gaussian()) + set_rescor(FALSE), data = dl)
writeLines(sc, file.path(out, "brms-mvln-stancode.stan"))
lines <- strsplit(sc, "\n")[[1]]
cat("lines with Llncor or nlncor or lncor:\n")
cat(paste0(grep("lncor", lines), ": ", grep("lncor", lines, value = TRUE)), sep = "\n")

## --- fits ------------------------------------------------------------
run <- function(nm, ...) {
  sec(nm)
  t0 <- Sys.time()
  fit <- tryCatch(do.call(brm, c(list(...), ctl)),
                  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(NULL)
  cat("seconds", round(as.numeric(Sys.time() - t0, units = "secs")), "\n")
  saveRDS(fit, file.path(out, paste0("brms-", nm, ".rds")))
  cat("max Rhat", max(rhat(fit), na.rm = TRUE), "\n")
  fit
}
post_cat <- function(fit, newg) {
  v <- variables(fit)
  cat("variables (b_, sd_, cor_):\n")
  print(grep("^(b_|sd_|cor_)", v, value = TRUE))
  print(posterior_summary(fit, variable = grep("^(b_|sd_|cor_)", v, value = TRUE)))
  try_shape("ranef", ranef(fit)$g)
  try_shape("coef", coef(fit)$g)
  vc <- VarCorr(fit)$g
  cat("VarCorr names:", names(vc), "\n")
  try_shape("VarCorr$g$sd", vc$sd)
  try_shape("VarCorr$g$cor", vc$cor)
  try_shape("fitted", fitted(fit))
  try_shape("fitted linear", fitted(fit, scale = "linear"))
  try_shape("predict", predict(fit))
  try_shape("posterior_linpred", posterior_linpred(fit))
  try_shape("posterior_epred", posterior_epred(fit))
  try_shape("posterior_predict", posterior_predict(fit))
  try_shape("log_lik", log_lik(fit))
  try_shape("fitted newdata new level",
            fitted(fit, newdata = newg, allow_new_levels = TRUE))
  try_shape("fitted newdata new level re_formula NA",
            fitted(fit, newdata = newg, re_formula = NA))
  print(round(fitted(fit, newdata = newg, re_formula = NA)[, "Estimate", ], 5))
  ce <- tryCatch(conditional_effects(fit, "x", categorical = TRUE)[[1]],
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(ce)) cat("ce:", ce, "\n") else {
    cat("ce columns:", names(ce), "\n"); cat("ce cats__ levels:", levels(ce$cats__), "\n")
    cat("ce dim:", dim(ce), "\n")
  }
}
d3 <- rev_cat_data(21, 3)
d4 <- rev_cat_data(22, 4)
newg <- data.frame(x = c(-1, 0, 1), g = c("g1", "g2", "new"))
f3 <- run("cat3-refb", bf(y ~ x + (1 | p | g)), data = d3,
          family = categorical(refcat = "b"))
if (!is.null(f3)) post_cat(f3, newg)
f4 <- run("cat4-refb", bf(y ~ x + (1 | p | g)), data = d4,
          family = categorical(refcat = "b"))
if (!is.null(f4)) post_cat(f4, newg)

fm <- run("mult-vt", bf(Yc | trials(N) ~ x), data = dm, family = multinomial())
if (!is.null(fm)) {
  print(fixef(fm))
  ft <- fitted(fm)
  shape("fitted", ft)
  print(round(ft[1:6, , 1], 4)); print(round(ft[1:6, , 3], 4))
  cat("trials of rows 1..6:", dm$N[1:6], "\n")
  cat("rows with N == 1:", head(which(dm$N == 1)), "\n")
  i1 <- which(dm$N == 1)[1:3]
  print(round(ft[i1, , 3], 5))
  try_shape("fitted linear", fitted(fm, scale = "linear"))
  try_shape("posterior_linpred", posterior_linpred(fm))
}

## --- lane's saved fits: refcat = NA and the mixture -------------------
sec("lane rds: refcat NA, mixture")
for (nm in c("cat-refna", "mix")) {
  p <- file.path("dev/catmulti-out", paste0("brms-", nm, ".rds"))
  if (!file.exists(p)) { cat(nm, "missing\n"); next }
  f <- readRDS(p)
  cat(nm, "family:", format(f$family$family), "\n")
  try_shape(paste(nm, "posterior_linpred"), posterior_linpred(f))
  try_shape(paste(nm, "fitted linear"), fitted(f, scale = "linear"))
  try_shape(paste(nm, "linpred transform"), posterior_linpred(f, transform = TRUE))
}
cat("DONE\n")
