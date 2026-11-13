# Reviewer of lane priortab: brms 2.23.0 default_prior() against frmtmb's
# default_prior() (route fit and sample), row by row on the slot key,
# lb and ub compared. Cases the lane's dev/priortab-brms.R does not cover.
# Rscript dev/reviews/priortab-review-compare.R > dev/reviews/priortab-review-compare.txt
lib <- Sys.getenv("PRIORTAB_LIB", "C:/Users/adf44/source/r/wt-priortab-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  loadNamespace("brms"); loadNamespace("frmtmb")
  loadNamespace("frmtmb.sample"); loadNamespace("frmtmb.spatial")
})
cat("frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("brms", as.character(packageVersion("brms")), "\n")
options(width = 160, warn = 1)
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$y2 <- d$y + rnorm(n)
d$yp <- abs(d$y) + 0.1
d$cnt <- rpois(n, 3)
d$pb <- plogis(d$y)
d$o <- cut(d$y, 4, labels = FALSE)
d$om <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d$xm <- d$x + rnorm(n, sd = 0.1)
d$sdx <- 0.1
d$sei <- 0.2
d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
d$ang <- (d$pb - 0.5) * 3
d$zo <- pmin(pmax(d$pb + sample(c(-1, 0, 1), n, TRUE, c(.1, .8, .1)), 0), 1)
d$tr <- 10
d$cb <- rbinom(n, 10, 0.4)
W <- matrix(0, 10, 10)
for (i in 1:9) W[i, i + 1] <- W[i + 1, i] <- 1
rownames(W) <- colnames(W) <- levels(d$g)

cases <- list(
  mv_rescor_dpar = quote(default_prior(
    mvbf(bf(y ~ z + (1 | g), sigma ~ x), bf(y2 ~ z + (1 | g)),
         rescor = TRUE), data = d)),
  mv_student_rescor = quote(default_prior(
    mvbf(bf(y ~ z), bf(y2 ~ z), rescor = TRUE), data = d,
    family = student())),
  mv_pcor = quote(default_prior(
    mvbf(bf(y ~ z + (1 + z | p | g)), bf(y2 ~ z + (1 | p | g))), data = d)),
  nl = quote(default_prior(
    bf(yp ~ a * exp(b * z), a ~ 1 + (1 | g), b ~ 1 + x, nl = TRUE),
    data = d)),
  nl_sigma = quote(default_prior(
    bf(yp ~ a * exp(b * z), a ~ 1, b ~ 1, sigma ~ x, nl = TRUE),
    data = d)),
  sigma_x = quote(default_prior(bf(y ~ z, sigma ~ x + (1 | g)), data = d)),
  sigma_identity = quote(default_prior(y ~ z, data = d,
                                       family = gaussian(link_sigma = "identity"))),
  student_nu_identity = quote(default_prior(
    y ~ z, data = d, family = student(link_nu = "identity"))),
  mix2 = quote(default_prior(y ~ z, data = d,
                             family = mixture(gaussian, gaussian))),
  mix_theta = quote(default_prior(bf(y ~ z, theta2 ~ x), data = d,
                                  family = mixture(gaussian, gaussian))),
  mix_pois_gamma = quote(default_prior(cnt ~ z, data = d,
                                       family = mixture(poisson, negbinomial))),
  cumulative = quote(default_prior(o ~ z + (1 | g), data = d,
                                   family = cumulative())),
  cumulative_disc = quote(default_prior(bf(o ~ z, disc ~ x), data = d,
                                        family = cumulative())),
  acat_cs = quote(default_prior(o ~ cs(z), data = d, family = acat())),
  smooth = quote(default_prior(y ~ s(x) + (1 | g), data = d)),
  smooth_sigma = quote(default_prior(bf(y ~ s(x), sigma ~ s(z)), data = d)),
  gp = quote(default_prior(y ~ gp(x), data = d)),
  gp_by = quote(default_prior(y ~ gp(x, by = cat), data = d)),
  mo = quote(default_prior(y ~ mo(om), data = d)),
  me = quote(default_prior(y ~ me(xm, sdx), data = d)),
  cor_slopes = quote(default_prior(y ~ z + (1 + z | g) + (1 | cat), data = d)),
  cor_dpar = quote(default_prior(bf(y ~ (1 + z | g), sigma ~ (1 + z | g)),
                                 data = d)),
  ar1 = quote(default_prior(y ~ z + ar(p = 1, gr = g), data = d)),
  arma = quote(default_prior(y ~ z + arma(p = 1, q = 1, gr = g), data = d)),
  ar_cov_pois = quote(default_prior(cnt ~ z + ar(p = 1, gr = g, cov = TRUE),
                                    data = d, family = poisson())),
  ma2 = quote(default_prior(y ~ z + ma(q = 2, gr = g), data = d)),
  cosy = quote(default_prior(y ~ z + cosy(gr = g), data = d)),
  se_sigma = quote(default_prior(bf(y | se(sei, sigma = TRUE) ~ z), data = d)),
  car = quote(default_prior(y ~ z + car(W, gr = g), data = d,
                            data2 = list(W = W))),
  weibull = quote(default_prior(yp ~ z, data = d, family = weibull())),
  exgauss = quote(default_prior(y ~ z, data = d, family = exgaussian())),
  vonmises = quote(default_prior(ang ~ z, data = d, family = von_mises())),
  gev = quote(default_prior(y ~ z, data = d, family = gen_extreme_value())),
  frechet = quote(default_prior(yp ~ z, data = d, family = frechet())),
  asym_laplace = quote(default_prior(y ~ z, data = d, family = asym_laplace())),
  skew_normal = quote(default_prior(y ~ z, data = d, family = skew_normal())),
  beta_binomial = quote(default_prior(cb | trials(tr) ~ z, data = d,
                                      family = beta_binomial())),
  hurdle_gamma = quote(default_prior(yp ~ z, data = d, family = hurdle_gamma())),
  zoib = quote(default_prior(zo ~ z, data = d,
                             family = zero_one_inflated_beta())),
  zib = quote(default_prior(zo ~ z, data = d, family = zero_inflated_beta())),
  categorical = quote(default_prior(cat ~ z, data = d, family = categorical())),
  shifted_ln = quote(default_prior(yp ~ z, data = d,
                                   family = shifted_lognormal())),
  shifted_ln_ndt = quote(default_prior(bf(yp ~ z, ndt ~ x), data = d,
                                       family = shifted_lognormal())),
  zi_x = quote(default_prior(bf(cnt ~ z, zi ~ x), data = d,
                             family = zero_inflated_poisson())),
  hurdle_pois = quote(default_prior(cnt ~ z, data = d,
                                    family = hurdle_poisson())),
  custom = quote(default_prior(cnt ~ z, data = d, family = custom_family(
    "myfam", dpars = c("mu", "phi"), links = c("log", "log"),
    lb = c(NA, 0), type = "int")))
)

norm_b <- function(v) {
  v <- as.character(v)
  v[is.na(v) | v == ""] <- "."
  num <- suppressWarnings(as.numeric(v))
  ok <- !is.na(num)
  v[ok] <- format(num[ok], digits = 7)
  trimws(v)
}
keycols <- c("class", "coef", "group", "resp", "dpar", "nlpar")
tidy <- function(tab) {
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  for (k in keycols) {
    if (!k %in% names(tab)) tab[[k]] <- ""
    tab[[k]][is.na(tab[[k]])] <- ""
  }
  tab$key <- do.call(paste, c(tab[keycols], sep = "|"))
  tab$lb <- norm_b(tab$lb)
  tab$ub <- norm_b(tab$ub)
  tab$prior <- as.character(tab$prior)
  tab$prior[tab$prior %in% c("", "(flat)")] <- "(flat)"
  tab
}
run <- function(expr, ns, extra = list()) {
  env <- list2env(c(list(d = d, W = W), extra), parent = asNamespace(ns))
  w <- character(0)
  r <- tryCatch(withCallingHandlers(eval(expr, env), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
    error = function(e) structure(conditionMessage(e), class = "err"))
  attr(r, "warns") <- w
  r
}
compare <- function(nm, b, f, route) {
  if (inherits(f, "err")) {
    cat(sprintf("  [%s] frmtmb ERROR: %s\n", route, substr(f, 1, 200)))
    return(invisible())
  }
  f <- tidy(f)
  only_b <- setdiff(b$key, f$key)
  only_f <- setdiff(f$key, b$key)
  for (k in only_b) cat(sprintf("  [%s] only brms : %-40s lb=%s ub=%s %s\n",
                                 route, k, b$lb[b$key == k][1],
                                 b$ub[b$key == k][1], b$prior[b$key == k][1]))
  for (k in only_f) cat(sprintf("  [%s] only frm  : %-40s lb=%s ub=%s %s\n",
                                 route, k, f$lb[f$key == k][1],
                                 f$ub[f$key == k][1], f$prior[f$key == k][1]))
  com <- intersect(b$key, f$key)
  for (k in com) {
    bi <- which(b$key == k)[1]; fi <- which(f$key == k)[1]
    if (b$lb[bi] != f$lb[fi] || b$ub[bi] != f$ub[fi]) {
      cat(sprintf("  [%s] BOUND DIFF: %-40s brms lb=%s ub=%s | frm lb=%s ub=%s\n",
                  route, k, b$lb[bi], b$ub[bi], f$lb[fi], f$ub[fi]))
    }
    if (route == "sample" && b$prior[bi] != f$prior[fi]) {
      cat(sprintf("  [%s] prior diff: %-40s brms %s | frm %s\n", route, k,
                  b$prior[bi], f$prior[fi]))
    }
  }
}
for (nm in names(cases)) {
  cat("\n=== ", nm, ": ", deparse1(cases[[nm]]), "\n", sep = "")
  b <- run(cases[[nm]], "brms")
  if (inherits(b, "err")) {
    cat("  brms ERROR:", substr(b, 1, 200), "\n"); next
  }
  if (length(attr(b, "warns"))) cat("  brms WARN:", attr(b, "warns"), "\n")
  b <- tidy(b)
  for (rt in c("fit", "sample")) {
    ex <- cases[[nm]]
    ex$route <- rt
    f <- run(ex, "frmtmb")
    compare(nm, b, f, rt)
  }
}
cat("\ndone\n")
