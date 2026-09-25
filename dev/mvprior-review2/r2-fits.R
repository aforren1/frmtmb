# Real fits on reviewer data. Per arm: fits with accepted priors, saved
# for identical() between arms; the sd-per-dpar targeting check on a
# shared (1 | ID | g) block; sample-route defaults.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-prelude.R")
suppressPackageStartupMessages({library(frmtmb.sample); library(frmtmb.latent)})
source(file.path(R2_ROOT, "dev/mvprior-review2/r2-cases.R"))
d <- r2_data()
out <- list()
fitsave <- function(nm, expr) {
  f <- tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) e)
  if (inherits(f, "error")) { cat(nm, "ERROR", conditionMessage(f), "\n"); out[[nm]] <<- conditionMessage(f); return(invisible(NULL)) }
  v <- list(ll = logLik(f), fixef = fixef(f), vc = tryCatch(VarCorr(f), error = function(e) NULL),
            par = f$fit$par %||% f$opt$par, ps = tryCatch(as.data.frame(prior_summary(f)), error = function(e) NULL))
  out[[nm]] <<- v
  cat(nm, "logLik", format(as.numeric(v$ll), digits = 12), "\n")
  invisible(f)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
cat4 <- categorical()
fitsave("cat4_b_mub", frm(bf(cat4 ~ x + z), family = cat4, data = d,
  prior = set_prior("normal(0, 0.1)", class = "b", dpar = "mub") +
    set_prior("normal(0, 1)", class = "Intercept", dpar = "mud")))
fitsave("cat2_b", frm(bf(cat2 ~ x), family = categorical(), data = d,
  prior = set_prior("normal(0, 0.2)", class = "b", dpar = "muyes")))
fitsave("mix3", frm(bf(ym3 ~ x), family = mixture(gaussian(), gaussian(), gaussian()), data = d,
  prior = set_prior("normal(0, 1)", class = "Intercept", dpar = "mu2") +
    set_prior("student_t(3, 0, 2.5)", class = "sigma3") +
    set_prior("normal(0, 0.5)", class = "b", dpar = "mu1", coef = "x")))
fitsave("mixsig", frm(bf(ym ~ x, sigma2 ~ z, theta2 ~ w), family = mixture(gaussian(), gaussian()), data = d,
  prior = set_prior("normal(0, 1)", class = "b", dpar = "theta2") +
    set_prior("student_t(3, 0, 2.5)", class = "sigma1")))
fitsave("mvcat", frm(bf(yg ~ x, family = gaussian()) + bf(cat4 ~ x, family = categorical()) + set_rescor(FALSE),
  data = d, prior = set_prior("normal(0, 0.3)", class = "b", dpar = "mub", resp = "cat4") +
    set_prior("normal(0, 1)", class = "b", resp = "yg")))
fitsave("zip", frm(bf(cnt ~ x + (1 | g), zi ~ z), family = zero_inflated_poisson(), data = d,
  prior = set_prior("normal(0, 1)", class = "b") + set_prior("normal(0, 1)", class = "b", dpar = "zi") +
    set_prior("student_t(3, 0, 1)", class = "sd")))
fitsave("cumul", frm(bf(ord ~ x + z), family = cumulative(), data = d,
  prior = set_prior("normal(0, 0.5)", class = "b") + set_prior("normal(0, 3)", class = "Intercept")))
fitsave("sratio", frm(bf(ord ~ x), family = sratio(), data = d,
  prior = set_prior("normal(0, 0.5)", class = "b")))
fitsave("hln", frm(bf(pos ~ x), family = hurdle_lognormal(), data = d,
  prior = set_prior("normal(0, 0.5)", class = "b") + set_prior("beta(2, 2)", class = "hu")))
fitsave("sigslope", frm(bf(yg ~ 1, sigma ~ x), data = d,
  prior = set_prior("normal(0, 0.2)", class = "b", dpar = "sigma")))

# (1 | ID | g) across mub, muc, mud: which SD does sd dpar = muc reach?
cat("\n-- sd targeting on (1 | ID | g) --\n")
for (dp in c("none", "mub", "muc", "mud")) {
  pr <- if (dp == "none") NULL else set_prior("normal(0, 0.02)", class = "sd", dpar = dp)
  f <- tryCatch(suppressWarnings(suppressMessages(frm(bf(cat4 ~ x + (1 | ID | g)), family = categorical(),
                                                      data = d, prior = pr))), error = function(e) e)
  if (inherits(f, "error")) { cat(dp, "ERROR", conditionMessage(f), "\n"); next }
  vc <- VarCorr(f)
  sds <- unlist(lapply(vc, function(v) if (is.list(v)) v$sd[, 1] else attr(v, "stddev")))
  cat(sprintf("%-5s", dp), paste(names(sds), format(sds, digits = 4), collapse = "  "), "\n")
  out[[paste0("idsd_", dp)]] <- list(ll = logLik(f), sds = sds)
}
# lca with dpar-qualified prior (both arms accept)
set.seed(3)
n <- 200; cl <- rbinom(n, 1, 0.4) + 1
pr <- rbind(c(0.85, 0.80, 0.75, 0.90), c(0.15, 0.20, 0.25, 0.10))
Y <- matrix(0L, n, 4); for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
dd <- data.frame(x = rnorm(n)); dd$Y <- Y
fitsave("lca", frm(bf(Y ~ x), family = lca(K = 2), data = dd,
  prior = set_prior("normal(0, 0.5)", class = "b", dpar = "theta1")))

cat("\n-- sample-route default priors --\n")
show_def <- function(nm, f, fam, dat) {
  tab <- tryCatch(as.data.frame(default_prior(f, data = dat, family = fam, route = "sample")),
                  error = function(e) conditionMessage(e))
  cat("==", nm, "\n")
  if (is.data.frame(tab)) {
    tab <- tab[tab$prior != "(flat)" & tab$class != "theta", ]
    print(tab[, c("prior", "class", "coef", "group", "resp", "dpar")], row.names = FALSE)
  } else cat(tab, "\n")
}
show_def("cat4re", bf(cat4 ~ x + (1 | g)), categorical(), d)
show_def("mix3", bf(ym3 ~ x), mixture(gaussian(), gaussian(), gaussian()), d)
show_def("mixdiff_re", bf(ym ~ x + (1 | g)), mixture(gaussian(), student()), d)
show_def("lca", bf(Y ~ x), lca(K = 3), dd)
saveRDS(out, file.path(R2_ROOT, "dev/mvprior-review2", paste0("r2-fits-", r2_arm, ".rds")))
