# Lane meme, fix round 1: short repros of the four pre-existing items
# the review filed (dev/reviews/2026-10-28-meme.md, "Pre-existing").
# Output: dev/meme-out/filed.txt
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 40
d <- data.frame(x = rnorm(n), z = rnorm(n), sx = 0.3, sz = 0.3,
                g = factor(rep(1:8, each = 5)))
d$xg <- rnorm(8)[d$g]; d$wg <- rnorm(8)[d$g]; d$sg <- 0.3
d$y <- 1 + d$x + rnorm(n)

cat("== 1. update(<draws>, prior = ) replaces the first call's prior\n")
f <- bf(y ~ me(x, sx) + me(z, sz))
ds <- q(frm_sample(f, data = d, family = gaussian(), chains = 1,
                   iter = 200, seed = 1, refresh = 0,
                   prior = set_prior("normal(0, 1)", class = "meanme") +
                     set_prior("exponential(1)", class = "sdme")))
ds2 <- q(update(ds, prior = set_prior("lkj(4)", class = "corme")))
ps <- as.data.frame(prior_summary(ds2))
print(ps[ps$class %in% c("meanme", "sdme", "corme"), c("prior", "class")])
cat("brms merges: meanme normal(0, 1) and sdme exponential(1) would stay\n")

cat("\n== 2. variables() orders me() hyperparameters by term\n")
fit <- frm(y ~ me(xg, sg, gr = g) + me(z, sz) + me(wg, sg, gr = g), data = d)
print(grep("^meanme", variables(fit), value = TRUE))
cat("brms: meanme_mexg meanme_mewg meanme_mez (by group)\n")

cat("\n== 3. default_prior() gives sigma no lb = 0\n")
tab <- as.data.frame(default_prior(y ~ x, data = d))
print(tab[tab$class == "sigma", c("prior", "class", "lb", "ub")])
cat("brms: sigma lb = 0\n")

cat("\n== 4. mi() latent values are miss_<k> in draws\n")
dm <- d
dm$x[c(2, 7)] <- NA
dsm <- q(frm_sample(bf(y ~ mi(x)) + bf(x | mi() ~ 1), data = dm,
                    family = gaussian(), chains = 1, iter = 200, seed = 1,
                    refresh = 0))
print(grep("miss|Ymi", colnames(dsm$draws), value = TRUE))
cat("brms: Ymi_x[1], Ymi_x[2]\n")
