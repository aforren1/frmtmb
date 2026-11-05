# brms 2.23.0's rules for class "sds": what stancode() accepts and writes.
# Seed 43. brms needs no compile for stancode().
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(43)
n <- 100
dd <- data.frame(x = runif(n), z = runif(n),
                 f = factor(rep(c("a", "b"), n / 2)),
                 g = factor(rep(1:10, each = 10)),
                 w = rnorm(n))
dd$y <- rnorm(n)
dd$y2 <- rnorm(n)
try_code <- function(lab, fo, pr, ...) {
  cat("\n---", lab, "\n")
  code <- tryCatch(suppressWarnings(brms::stancode(fo, data = dd, prior = pr,
                                                   ...)),
                   error = function(e) e)
  if (inherits(code, "error")) {
    cat("ERROR:", substr(gsub("[\r\n]+", " ", conditionMessage(code)), 1,
                         300), "\n")
    return(invisible())
  }
  ln <- strsplit(as.character(code), "\n")[[1]]
  cat(grep("sds_|sd_[0-9]+.*target|real<lower=0> sds|vector<lower=0>",
           ln, value = TRUE), sep = "\n")
  invisible()
}
try_code("sds class", y ~ s(x) + s(z), prior(normal(0, 1), class = sds))
try_code("sds coef exact", y ~ s(x) + s(z),
         prior(normal(0, 1), class = sds, coef = "s(x)"))
try_code("sds coef both", y ~ s(x) + s(z),
         prior(normal(0, 1), class = sds) +
           prior(cauchy(0, 2), class = sds, coef = "s(z)"))
try_code("sds coef no space", y ~ s(z, k = 5),
         prior(normal(0, 1), class = sds, coef = "s(z,k=5)"))
try_code("sds coef with space", y ~ s(z, k = 5),
         prior(normal(0, 1), class = sds, coef = "s(z, k = 5)"))
try_code("written tight", y ~ s(z,k=5),
         prior(normal(0, 1), class = sds, coef = "s(z, k = 5)"))
try_code("sds group", y ~ s(x),
         prior(normal(0, 1), class = sds, group = "s(x)"))
try_code("sd no smooth", y ~ s(x), prior(normal(0, 1), class = sd))
try_code("sd group smooth", y ~ s(x) + (1 | g),
         prior(normal(0, 1), class = sd, group = "s(x)"))
try_code("sd beside RE: s not touched", y ~ s(x) + (1 | g),
         prior(normal(0, 1), class = sd))
try_code("t2 coef", y ~ t2(x, z),
         prior(normal(0, 1), class = sds, coef = "t2(x, z)"))
try_code("by coef", y ~ f + s(x, by = f),
         prior(normal(0, 1), class = sds, coef = "s(x, by = f)"))
try_code("by level coef", y ~ f + s(x, by = f),
         prior(normal(0, 1), class = sds, coef = "s(x):fa"))
try_code("no smooth", y ~ x, prior(normal(0, 1), class = sds))
try_code("sds lb on class", y ~ s(x),
         prior(normal(0, 1), class = sds, lb = 0.1))
try_code("sds lb on coef", y ~ s(x),
         prior(normal(0, 1), class = sds, coef = "s(x)", lb = 0.1))
try_code("sds lb 0 on coef", y ~ s(x),
         prior(normal(0, 1), class = sds, coef = "s(x)", lb = 0))
try_code("dpar", bf(y ~ s(x), sigma ~ s(z)),
         prior(normal(0, 1), class = sds, dpar = "sigma"))
try_code("no dpar reaches mu only", bf(y ~ s(x), sigma ~ s(z)),
         prior(normal(0, 1), class = sds))
try_code("nlpar", bf(y ~ a * exp(b * x), a ~ s(z), b ~ 1, nl = TRUE),
         prior(normal(0, 1), class = sds, nlpar = "a") +
           prior(normal(0, 1), nlpar = "b"))
try_code("mv resp", bf(y ~ s(x)) + bf(y2 ~ s(z)) + set_rescor(FALSE),
         prior(normal(0, 1), class = sds, resp = "y"))
try_code("mv no resp", bf(y ~ s(x)) + bf(y2 ~ s(z)) + set_rescor(FALSE),
         prior(normal(0, 1), class = sds))
try_code("lkj", y ~ s(x), prior(lkj(1), class = sds))
cat("\n=== tables\n")
print(as.data.frame(default_prior(y ~ s(z, bs = "cr", k = 5) + s(x, z),
                                  data = dd))[, c("prior", "class", "coef",
                                                  "lb")])
print(as.data.frame(default_prior(bf(y ~ s(x)) + bf(y2 ~ s(z)) +
                                    set_rescor(FALSE), data = dd))[
  , c("prior", "class", "coef", "resp", "lb")])
print(as.data.frame(default_prior(bf(y ~ a * exp(b * x), a ~ s(z), b ~ 1,
                                     nl = TRUE), data = dd))[
  , c("prior", "class", "coef", "nlpar", "lb")])
print(as.data.frame(default_prior(y ~ s(x, w) + s(g, bs = "re") +
                                    s(x, g, bs = "fs", k = 4), data = dd))[
  , c("prior", "class", "coef", "lb")])
