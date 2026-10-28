# Lane reslot: what brms 2.23.0 means by a (flat) slot, run rather than
# recalled: standata() and stancode() only, no sampling. Seed 1 data.
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
G <- 8L; n <- 10L
dd <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n),
                 z = rnorm(G * n), t = rep(seq_len(n), G))
dd$y <- 1 + 0.5 * dd$x + rnorm(G)[dd$g] + rnorm(G * n)
dd$m <- rep(1:4, length.out = nrow(dd))
dd$ord <- cut(dd$y, 4, labels = FALSE)
grab <- function(code, pat) {
  l <- strsplit(as.character(code), "\n")[[1]]
  cat(paste0("  | ", grep(pat, l, value = TRUE)), sep = "\n")
}

cat("\n-- the delta default of an equidistant cumulative():\n")
print(default_prior(bf(ord ~ x), data = dd,
                    family = cumulative(threshold = "equidistant")))
grab(stancode(bf(ord ~ x), data = dd,
              family = cumulative(threshold = "equidistant")),
     "lprior|delta;")

cat("\n-- con_simo under simo \"\" and under the default:\n")
print(standata(bf(y ~ mo(m)), data = dd,
               prior = set_prior("", class = "simo", coef = "mom1"))$con_simo_1)
print(standata(bf(y ~ mo(m)), data = dd)$con_simo_1)

cat("\n-- con_theta of a mixture under theta \"\" and under the default:\n")
mp <- c(set_prior("normal(0, 5)", class = "Intercept", dpar = "mu1"),
        set_prior("normal(1, 5)", class = "Intercept", dpar = "mu2"))
print(standata(bf(y ~ x), data = dd, family = mixture(gaussian, gaussian),
               prior = c(mp, set_prior("", class = "theta")))$con_theta)
print(standata(bf(y ~ x), data = dd, family = mixture(gaussian, gaussian),
               prior = mp)$con_theta)

cat("\n-- lb on class cor: does the Stan program carry it?\n")
grab(stancode(bf(y ~ x + (1 + x | g)), data = dd,
              prior = set_prior("lkj(2)", class = "cor", lb = 0)),
     "cholesky_factor_corr|lkj")
grab(stancode(bf(y ~ x + (1 + x | g)), data = dd,
              prior = set_prior("", class = "cor", group = "g", lb = -0.5,
                                ub = 0.5)),
     "cholesky_factor_corr|lkj")

cat("\n-- every slot (flat): which parameters keep a density\n")
fl <- c(set_prior("", class = "Intercept"), set_prior("", class = "sd"),
        set_prior("", class = "sigma"), set_prior("", class = "cor"),
        set_prior("", class = "ar"))
grab(stancode(bf(y ~ x + (1 + x + z | g) + ar(time = t, gr = g, cov = TRUE)),
              data = dd, prior = fl),
     "lprior|<lower|cholesky|ar;")
cat("\n-- the sigma brms's (flat) is on under ar(cov = TRUE):\n")
grab(stancode(bf(y ~ x + ar(time = t, gr = g, cov = TRUE)), data = dd,
              prior = c(set_prior("", class = "sigma"),
                        set_prior("", class = "Intercept"))),
     "normal_time_het|sigma")
