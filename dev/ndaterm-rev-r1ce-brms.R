# Reviewer of lane ndaterm, fix round 1: the values brms 2.23.0's
# conditional_effects() grid holds a vint() and a thres(gr = ) variable
# at (brms:::prepare_conditions() on empty fits, no compile). Seed 7.
# Log: dev/ndaterm-rev-out/r1ce-brms.txt (redirected by the caller)
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(7)
n <- 200
d <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE),
                gnum = sample(1:3, n, TRUE), up = rbinom(n, 1, 0.6),
                rt = 0.3 + rexp(n))
d$y <- rbinom(n, d$N, plogis(0.3 + 0.5 * d$x))
d$yo <- sample(1:4, n, TRUE)
d$gf <- factor(d$gnum)
cbin <- custom_family("cbin", dpars = "mu", links = "logit", type = "int",
                      vars = "vint1[n]")
fv <- brm(bf(y | vint(N) ~ x), family = cbin, data = d, empty = TRUE)
cv <- brms:::prepare_conditions(fv, effects = "x")
cat("vint(N): grid N =", cv$N, "; median(N) =", median(d$N),
    "mean(N) =", mean(d$N), "\n")
ft <- brm(bf(yo | thres(gr = gf) ~ x), family = cumulative(), data = d,
          empty = TRUE)
ct <- brms:::prepare_conditions(ft, effects = "x")
cat("thres(gr = gf): grid gf =", as.character(ct$gf), "; levels", levels(d$gf),
    "\n")
fw <- brm(bf(rt | dec(up) ~ x), family = wiener(), data = d, empty = TRUE)
cw <- brms:::prepare_conditions(fw, effects = "x")
cat("dec(up): grid up =", cw$up, "; mean =", mean(d$up), "\n")
fb <- brm(bf(y | trials(N) ~ x), family = binomial(), data = d,
          empty = TRUE)
cb <- brms:::prepare_conditions(fb, effects = "x")
cat("trials(N): grid N =", cb$N, "\n")
fbx <- brm(bf(y | vint(N) ~ x + N), family = cbin, data = d, empty = TRUE)
cbx <- brms:::prepare_conditions(fbx, effects = "x")
cat("vint(N) with N also a predictor: grid N =", cbx$N, "\n")
cat("done\n")
