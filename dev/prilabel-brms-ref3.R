# Lane prilabel: brms 2.23.0's per-coefficient sd rows on gr(by = ),
# mo(), mmc(), a nonlinear parameter and an |ID| block spanning mu and
# sigma (no compile). Data: dev/dblbar-data.R, seed 20261013.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(brms))
dd <- dblbar_data()
set.seed(7)
dd$x1 <- rnorm(nrow(dd))
dd$x2 <- rnorm(nrow(dd))
sink(file.path(wt, "dev/prilabel-log/brms-ref3.txt"), split = TRUE)
cat("brms", format(packageVersion("brms")), "\n")
rows <- function(f, ...) {
  cat("\n== default_prior:", deparse1(f$formula), "\n")
  p <- as.data.frame(default_prior(f, data = dd, ...))
  print(p[p$class %in% c("sd"), c("class", "coef", "group", "dpar",
                                   "nlpar")], row.names = FALSE)
}
sdlines <- function(f, prior, ...) {
  r <- tryCatch({
    sc <- stancode(f, data = dd, prior = prior, ...)
    grep("lprior.*sd_|sd_[0-9]+\\[", strsplit(sc, "\n")[[1]], value = TRUE)
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(r, sep = "\n")
}
f1 <- bf(y ~ x + (1 + x | gr(g, by = hg)))
rows(f1)
cat("-- coef x on gr(g, by = hg)\n")
sdlines(f1, set_prior("normal(0, 0.5)", class = "sd", group = "g",
                      coef = "x"))
f2 <- bf(y ~ x + mo(xo) + (1 + mo(xo) | g))
rows(f2)
cat("-- coef moxo\n")
sdlines(f2, set_prior("normal(0, 0.5)", class = "sd", group = "g",
                      coef = "moxo"))
f3 <- bf(y ~ x + (1 + mmc(x1, x2) | mm(m1, m2)))
rows(f3)
f4 <- bf(y ~ a * exp(b * x), a ~ 1 + (1 + z | g), b ~ 1, nl = TRUE)
rows(f4, prior = c(prior(normal(1, 1), nlpar = "a"),
                   prior(normal(0, 1), nlpar = "b")))
cat("-- coef z on nlpar a\n")
sdlines(f4, c(prior(normal(1, 1), nlpar = "a"),
              prior(normal(0, 1), nlpar = "b"),
              set_prior("normal(0, 0.5)", class = "sd", group = "g",
                        coef = "z", nlpar = "a")))
f5 <- bf(y ~ x + (1 + x | p | g), sigma ~ (1 | p | g))
rows(f5)
cat("-- coef Intercept, group g, no dpar, on a block spanning mu, sigma\n")
sdlines(f5, set_prior("normal(0, 0.5)", class = "sd", group = "g",
                      coef = "Intercept"))
cat("-- coef Intercept, group g, dpar sigma\n")
sdlines(f5, set_prior("normal(0, 0.5)", class = "sd", group = "g",
                      coef = "Intercept", dpar = "sigma"))
cat("-- coef x with lb/ub only (no density)\n")
sdlines(bf(y ~ x + (1 + x | g)),
        set_prior("", class = "sd", group = "g", coef = "x", ub = 3))
cat("-- coef x, class sd, group g, and a class-wide sd prior\n")
sdlines(bf(y ~ x + (1 + x | g)),
        c(set_prior("normal(0, 2)", class = "sd"),
          set_prior("normal(0, 0.5)", class = "sd", group = "g",
                    coef = "x")))
f6 <- bf(y ~ x + (1 + x | g) + (1 | g1), sigma ~ (1 + x | g))
rows(f6)
mv <- bf(y ~ x + (1 + x | g)) + bf(z ~ x + (1 | g)) + set_rescor(FALSE)
cat("\n== mvbf default_prior\n")
p <- as.data.frame(default_prior(mv, data = dd))
print(p[p$class %in% c("sd"), c("class", "coef", "group", "resp")],
      row.names = FALSE)
cat("-- mv coef x, group g, resp y\n")
sdlines(mv, set_prior("normal(0, 0.5)", class = "sd", group = "g",
                      coef = "x", resp = "y"))
sink()
