# Lane meme, fix round 1: what brms 2.23.0 prints on group rows (sd,
# cor, corme, grouped thresholds, delta) and how it reads a coef with a
# bound on b, sd, meanme and sdme. Output: dev/meme-out/fix-brms-grp.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(3)
d <- data.frame(x = rnorm(60), w = rnorm(60), g = factor(rep(1:6, 10)),
                th = factor(rep(c("a", "b"), 30)), sx = 0.3)
d$xg <- rnorm(6)[d$g]; d$zg <- rnorm(6)[d$g]; d$sg <- 0.3
d$y <- rnorm(60); d$yo <- sample(1:4, 60, TRUE)
show <- function(t) print(as.data.frame(t)[, c("prior", "class", "coef", "group", "lb", "ub")])
cat("== sd/cor/corme group rows\n")
show(default_prior(bf(y ~ x + (1 + x | g) + me(xg, sg, gr = g) + me(zg, sg, gr = g)), data = d))
cat("== grouped thresholds, flexible\n")
show(default_prior(bf(yo | thres(gr = th) ~ x), data = d, family = cumulative()))
cat("== grouped thresholds, equidistant\n")
show(default_prior(bf(yo | thres(gr = th) ~ x), data = d, family = cumulative(threshold = "equidistant")))
cat("== coef + bound\n")
f <- bf(y ~ x + w + (1 + x | g) + me(xg, sg, gr = g) + me(zg, sg, gr = g))
for (p in list(set_prior("normal(0, 1)", class = "b", coef = "w", ub = 0.1),
               set_prior("normal(0, 1)", class = "sd", group = "g", coef = "x", ub = 2),
               set_prior("normal(0, 1)", class = "sd", group = "g", ub = 2),
               set_prior("normal(0, 1)", class = "meanme", coef = "mexg", lb = 0),
               set_prior("exponential(1)", class = "sdme", coef = "mexg", ub = 0.5),
               set_prior("normal(0, 1)", class = "b", lb = -1),
               set_prior("exponential(1)", class = "sdme", ub = 3),
               set_prior("", class = "sdme", coef = "mexg", ub = 3))) {
  cat("--", p$class, p$coef, p$group, "lb", p$lb, "ub", p$ub, ":\n")
  r <- tryCatch({ sc <- stancode(f, data = d, prior = p)
    l <- strsplit(sc, "\n")[[1]]; grep("upper|lower=-1|lower=0>[[]Mme|meanme_1;|sdme_1;", l, value = TRUE) },
    error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(r, sep = "\n")
  r2 <- tryCatch({ set_prior(p$prior, class = p$class, coef = p$coef, group = p$group, lb = p$lb, ub = p$ub); "set_prior ok" },
    error = function(e) paste("set_prior ERROR:", conditionMessage(e)))
  cat(r2, "\n")
}
cat("== me(x) without sdx\n")
show(default_prior(bf(y ~ me(x)), data = d))
print(tryCatch(validate_prior(set_prior("normal(0,1)", class = "meanme"), bf(y ~ me(x)), data = d)[, 1:4], error = function(e) conditionMessage(e)))
