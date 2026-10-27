# Reviewer: brms 2.23.0 on every prior route the review attacks.
# default_prior(), validate_prior() and stancode() only (no fit).
# Usage: Rscript dev/review-acprior-brms.R
# Log: dev/review-acprior-log/brms.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 40
d <- data.frame(y = rnorm(n), y1 = rnorm(n), y2 = rnorm(n), x = rnorm(n),
                t = rep(1:10, 4), g = rep(1:4, each = 10),
                yp = rexp(n) + 0.5)
tab <- function(x) {
  print(as.data.frame(x)[, c("prior", "class", "coef", "group", "resp",
                             "dpar", "lb", "ub", "source")])
}
try_ <- function(expr) {
  tryCatch(expr, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}
decl <- function(fo, prior = NULL, pat, data = d, family = gaussian()) {
  sc <- try_(stancode(fo, data = data, prior = prior, family = family))
  if (is.character(sc)) {
    ln <- strsplit(sc, "\n")[[1]]
    cat(paste0("  | ", grep(pat, ln, value = TRUE)), sep = "\n")
  }
}

cat("\n===== mv: ar cov in y1, cosy in y2\n")
fmv <- bf(y1 ~ x + ar(t, g, cov = TRUE)) + bf(y2 ~ x + cosy(t, g)) +
  set_rescor(FALSE)
tab(default_prior(fmv, data = d))
cat("-- class ar with no resp in mv:\n")
try_(print(validate_prior(set_prior("normal(0.5, 0.2)", class = "ar"),
                          fmv, data = d)))
cat("-- class ar resp y1, lb = 0:\n")
decl(fmv, set_prior("", class = "ar", resp = "y1", lb = 0),
     "ar_y1;|cosy_y2;|lprior.*(ar|cosy)")
cat("-- class cosy resp y2 normal(0, 0.1):\n")
decl(fmv, set_prior("normal(0, 0.1)", class = "cosy", resp = "y2"),
     "ar_y1;|cosy_y2;|lprior.*(ar|cosy)")

cat("\n===== mv: ar cov in BOTH\n")
fmv2 <- bf(y1 ~ x + ar(t, g, cov = TRUE)) + bf(y2 ~ x + ar(t, g, cov = TRUE)) +
  set_rescor(FALSE)
tab(default_prior(fmv2, data = d))
cat("-- class ar, no resp, lb = 0.2 (stancode):\n")
decl(fmv2, set_prior("", class = "ar", lb = 0.2),
     "ar_y[12];|lprior.*ar")
cat("-- class ar, no resp, normal(0.5, 0.2) (stancode):\n")
decl(fmv2, set_prior("normal(0.5, 0.2)", class = "ar"),
     "ar_y[12];|lprior.*ar")

cat("\n===== mv: cosy in BOTH, class cosy no resp normal(0, 0.1)\n")
fmv3 <- bf(y1 ~ x + cosy(t, g)) + bf(y2 ~ x + cosy(t, g)) + set_rescor(FALSE)
decl(fmv3, set_prior("normal(0, 0.1)", class = "cosy"),
     "cosy_y[12];|lprior.*cosy")

cat("\n===== autocor on a dpar\n")
try_(tab(default_prior(bf(y ~ x, sigma ~ ar(t, g)), data = d)))

cat("\n===== ar(gr=) only time, arma(p = 2, q = 1)\n")
tab(default_prior(y ~ x + arma(t, g, p = 2, q = 1), data = d))
cat("-- ar(p = 2, cov = TRUE):\n")
try_(default_prior(y ~ x + ar(t, g, p = 2, cov = TRUE), data = d))
cat("-- ar() without gr:\n")
tab(default_prior(y ~ x + ar(tt), data = transform(d, tt = seq_len(n))))

cat("\n===== prior_string and cosy lb = -0.1 bounds only\n")
decl(y ~ x + cosy(t, g), prior_string("normal(0, 0.5)", class = "cosy"),
     "cosy;|lprior.*cosy")
decl(y ~ x + cosy(t, g), set_prior("", class = "cosy", lb = -0.1),
     "cosy;|lprior.*cosy")
cat("-- cosy ub = 1.5 (past brms's declared end):\n")
decl(y ~ x + cosy(t, g), set_prior("", class = "cosy", ub = 1.5),
     "cosy;|lprior.*cosy")
cat("-- ar cov lb = -1 (the table's own value):\n")
decl(y ~ x + ar(t, g, cov = TRUE), set_prior("", class = "ar", lb = -1),
     "ar;|lprior.*ar")

cat("\n===== unstr: cortime lkj(2); cortime lb\n")
decl(y ~ x + unstr(t, g), set_prior("lkj(2)", class = "cortime"),
     "cortime|lprior.*cortime")
try_(validate_prior(set_prior("lkj(2)", class = "cortime", lb = 0),
                    y ~ x + unstr(t, g), data = d))

cat("\n===== class theta on an RE / rescor model\n")
try_(validate_prior(set_prior("normal(0, 1)", class = "theta"),
                    y ~ x + (1 | g), data = d))
tab(default_prior(y ~ x + (1 + x | g), data = d))

cat("\n===== bounded natural defaults under a user prior\n")
cat("-- sigma, identity link, normal(1, 1):\n")
decl(bf(y ~ x, sigma ~ 1), set_prior("normal(1, 1)", class = "Intercept",
                                      dpar = "sigma"),
     "Intercept_sigma;|lprior.*sigma", family = brmsfamily("gaussian", link_sigma = "identity"))
decl(y ~ x, set_prior("normal(1, 1)", class = "sigma"),
     "sigma;|lprior.*sigma", family = brmsfamily("gaussian", link_sigma = "identity"))
cat("-- nu, student, normal(5, 2):\n")
decl(y ~ x, set_prior("normal(5, 2)", class = "nu"), "nu;|lprior.*nu",
     family = student())
cat("-- ndt, shifted_lognormal, normal(0.2, 0.1):\n")
decl(yp ~ x, set_prior("normal(0.2, 0.1)", class = "ndt"),
     "ndt;|lprior.*ndt|min_Y", family = shifted_lognormal())
cat("-- sd, normal(0, 1):\n")
decl(y ~ x + (1 | g), set_prior("normal(0, 1)", class = "sd"),
     "sd_1;|lprior.*sd")
cat("-- sigma lb = 0.5 bounds only:\n")
decl(y ~ x, set_prior("", class = "sigma", lb = 0.5), "sigma;|lprior.*sigma")
