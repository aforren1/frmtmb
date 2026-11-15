# Review of lane silent2: attacks on 5235, 5246, 5353 (core).
# Usage: Rscript dev/reviews/silent2-review-core.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-silent2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", find.package("frmtmb"), "\n")
try_s <- function(e) {
  ws <- character()
  r <- tryCatch(withCallingHandlers({force(e); "ran"},
                  warning = function(w) {
                    ws <<- c(ws, conditionMessage(w))
                    invokeRestart("muffleWarning")
                  }, message = function(m) invokeRestart("muffleMessage")),
                error = function(err) paste("ERROR:",
                  substr(gsub("\n", " ", conditionMessage(err)), 1, 160)))
  paste0(r, if (length(ws)) paste0(" [", length(ws), " warn: ",
                                   substr(ws[1], 1, 100), "]"))
}
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
say <- function(lab, e) cat(sprintf("%-60s %s\n", lab, try_s(e)))

## ---- 5246 anova: false refusals ----
set.seed(11)
n <- 80
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:8, 10)))
d$y <- d$x + rnorm(n)
d$cnt <- rpois(n, exp(0.5 + 0.3 * d$x))
d$nt <- 10L
d$s <- rbinom(n, 10, plogis(0.3 * d$x))
d$fl <- d$nt - d$s
d$w <- runif(n, 0.5, 2)
say("anova y ~ 1 vs I(y) ~ x", anova(q(y ~ 1, data = d), q(I(y) ~ x, data = d)))
say("anova poisson vs negbinomial (same y)",
    anova(q(cnt ~ x, family = poisson(), data = d),
          q(cnt ~ x, family = negbinomial(), data = d)))
say("anova poisson y~1 vs y~x",
    anova(q(cnt ~ 1, family = poisson(), data = d),
          q(cnt ~ x, family = poisson(), data = d)))
say("anova binomial s|trials(nt) ~1 vs ~x",
    anova(q(s | trials(nt) ~ 1, family = binomial(), data = d),
          q(s | trials(nt) ~ x, family = binomial(), data = d)))
say("anova binomial trials(nt) vs trials(10)",
    anova(q(s | trials(nt) ~ 1, family = binomial(), data = d),
          q(s | trials(10) ~ x, family = binomial(), data = d)))
say("anova weights vs no weights",
    anova(q(y ~ 1, data = d), q(y | weights(w) ~ x, data = d)))
say("anova gaussian vs student (same y)",
    anova(q(y ~ x, data = d), q(y ~ x, family = student(), data = d)))
say("anova mixed (1|g) vs none",
    anova(q(y ~ x, data = d), q(y ~ x + (1 | g), data = d)))
# NA rows dropped the same way in both, rows shuffled in one
dn <- d
dn$x[c(30, 31, 55)] <- NA
perm <- sample(n)
say("anova NA dropped by fit, one copy shuffled",
    anova(q(y ~ x, data = dn), q(y ~ x + z, data = dn[perm, ])))
say("anova NA dropped by fit vs complete.cases copy",
    anova(q(y ~ x, data = dn),
          q(y ~ x + z, data = dn[complete.cases(dn), ])))
dcc <- dn[complete.cases(dn), ]
rownames(dcc) <- NULL
say("anova NA dropped by fit vs complete.cases, names reset, shuffled",
    anova(q(y ~ x, data = dn), q(y ~ x + z, data = dcc[sample(nrow(dcc)), ])))
# character row names
dch <- d
rownames(dch) <- paste0("r", seq_len(n))
say("anova character row names, one shuffled",
    anova(q(y ~ 1, data = dch), q(y ~ x, data = dch[perm, ])))
# tibble / data.table
if (requireNamespace("tibble", quietly = TRUE)) {
  tb <- tibble::as_tibble(dn)
  say("anova tibble with NA, fit vs filtered tibble",
      anova(q(y ~ x, data = tb), q(y ~ x + z, data = tb[!is.na(tb$x), ])))
}
if (requireNamespace("data.table", quietly = TRUE)) {
  dt <- data.table::as.data.table(dn)
  say("anova data.table, NA dropped",
      anova(q(y ~ x, data = dt), q(y ~ x + z, data = dt)))
}
# subset() addition term vs pre-subsetted data
d$keep <- seq_len(n) > 10
say("anova subset(keep) vs data[keep, ]",
    anova(q(y | subset(keep) ~ 1, data = d), q(y ~ x, data = d[d$keep, ])))
# multivariate
say("anova mv rescor FALSE nested",
    anova(q(bf(y ~ 1) + bf(z ~ 1) + set_rescor(FALSE), data = d),
          q(bf(y ~ x) + bf(z ~ x) + set_rescor(FALSE), data = d)))
# ordinal: integer vs ordered factor response
d$o <- cut(d$y, c(-Inf, -1, 0, 1, Inf), labels = FALSE)
d$of <- factor(d$o, ordered = TRUE)
say("anova cumulative int y ~1 vs ~x",
    anova(q(o ~ 1, family = cumulative(), data = d),
          q(o ~ x, family = cumulative(), data = d)))
say("anova cumulative int vs ordered factor y",
    anova(q(o ~ 1, family = cumulative(), data = d),
          q(of ~ x, family = cumulative(), data = d)))
# bernoulli integer vs logical response
d$b <- as.integer(d$y > 0)
d$bl <- d$y > 0
say("anova bernoulli int vs logical y",
    anova(q(b ~ 1, family = bernoulli(), data = d),
          q(bl ~ x, family = bernoulli(), data = d)))
# large-scale gaussian (autoscale)
d$yb <- 1e6 * d$y + 5e7
say("anova gaussian y*1e6 ~1 vs ~x",
    anova(q(yb ~ 1, data = d), q(yb ~ x, data = d)))
# cens: same y, different censoring indicator (other likelihood data)
d$c1 <- ifelse(d$y > 1, "right", "none")
d$c2 <- "none"
say("anova cens(c1) vs cens(c2) (different data, should refuse?)",
    anova(q(y | cens(c1) ~ 1, data = d), q(y | cens(c2) ~ x, data = d)))
# mi(): response with NA modeled vs NA rows dropped
dm <- d
dm$y[1:5] <- NA
say("anova y | mi() ~ x vs y | mi() ~ x + z",
    anova(q(bf(y | mi() ~ x), data = dm), q(bf(y | mi() ~ x + z), data = dm)))

## ---- 5235: in-sample and constant payloads ----
scaled <- frmtmb_family(
  "scaled_normal", dpars = c("mu", "sigma"),
  links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms) {
    m <- dpars[["mu"]] * aterms[["vreal1"]]
    -0.5 * ((y - m) / dpars[["sigma"]])^2 - log(dpars[["sigma"]]) -
      0.5 * log(2 * pi)
  },
  init_dpars = list(mu = function(y, aterms) mean(y / aterms[["vreal1"]]),
                    sigma = function(y, aterms) stats::sd(y)),
  post = list(mean_fn = function(dpars, aterms) dpars[["mu"]] * aterms[["vreal1"]]),
  sim = function(dpars, aterms, n) {
    stats::rnorm(n, dpars[["mu"]] * aterms[["vreal1"]], dpars[["sigma"]])
  })
set.seed(3)
n <- 120
da <- data.frame(x = rnorm(n))
da$s <- rep_len(c(1, 3), n)
da$y <- (1 + 0.5 * da$x) * da$s + rnorm(n, 0, 0.3)
fe <- q(bf(y | vreal(rep_len(c(1, 3), 120)) ~ x), family = scaled, data = da)
say("5235 fitted()", fitted(fe))
say("5235 predict()", predict(fe))
say("5235 residuals()", residuals(fe))
say("5235 log_lik()", log_lik(fe))
say("5235 posterior_predict()", posterior_predict(fe, ndraws = 5))
say("5235 simulate()", simulate(fe, nsim = 2))
say("5235 fitted(newdata = data)", fitted(fe, newdata = da))
say("5235 loo()", loo(fe))
fk <- q(bf(y | vreal(rep(2, 120)) ~ x), family = scaled, data = da)
say("5235 constant rep(2,120) newdata 120", fitted(fk, newdata = da))
say("5235 constant rep(2,120) CE res 120",
    conditional_effects(fk, effects = "x", resolution = 120))
say("5235 constant rep(2,120) CE default", conditional_effects(fk, effects = "x"))
db <- data.frame(x = rnorm(60))
db$s <- rbinom(60, 10, plogis(db$x))
fb <- q(s | trials(10) ~ x, family = binomial(), data = db)
say("5235 trials(10) newdata", fitted(fb, newdata = db[1:7, ]))
say("5235 trials(10) CE", conditional_effects(fb, effects = "x"))
db$nn <- 10L
fb2 <- q(s | trials(nn) ~ x, family = binomial(), data = db)
say("5235 trials(nn) CE", conditional_effects(fb2, effects = "x"))
db$off <- runif(60, 1, 3)
db$cnt <- rpois(60, db$off * exp(0.2 * db$x))
fr <- q(cnt | rate(off * 2) ~ x, family = poisson(), data = db)
say("5235 rate(off * 2) newdata", fitted(fr, newdata = db[1:7, ]))
say("5235 rate(off * 2) CE", conditional_effects(fr, effects = "x"))
db$sey <- runif(60, 0.5, 1)
db$yy <- db$x + rnorm(60, 0, db$sey)
fse <- q(yy | se(sey) ~ x, data = db)
say("5235 se(sey) newdata", fitted(fse, newdata = db[1:7, ]))
say("5235 se(sey) CE predict", conditional_effects(fse, effects = "x", method = "predict"))

## ---- 5353 contrasts ----
set.seed(7)
dc <- data.frame(f = factor(sample(c("a", "b", "c"), 90, TRUE)),
                 o = factor(sample(c("lo", "mid", "hi"), 90, TRUE),
                            levels = c("lo", "mid", "hi"), ordered = TRUE),
                 g = factor(rep(1:10, each = 9)), x = rnorm(90),
                 ch = sample(c("p", "q"), 90, TRUE),
                 cyl = sample(c(4, 6, 8), 90, TRUE))
dc$g2 <- factor(rep(c(2:10, 1), each = 9), levels = 1:10)
dc$y <- as.numeric(dc$f) + rnorm(10)[dc$g] + rnorm(90)
nd <- dc[c(1:6, 30), ]
cmp <- function(lab, fml, set_opt, reset_opt, ...) {
  r1 <- withr::with_options(list(contrasts = set_opt), {
    f <- q(fml, data = dc, ...)
    list(f = f, p = tryCatch(fitted(f, newdata = nd)[, "Estimate"],
                             error = function(e) conditionMessage(e)),
         ins = fitted(f)[c(1:6, 30), "Estimate"])
  })
  p2 <- withr::with_options(list(contrasts = reset_opt),
          tryCatch(fitted(r1$f, newdata = nd)[, "Estimate"],
                   error = function(e) conditionMessage(e)))
  if (is.character(r1$p) || is.character(p2)) {
    cat(sprintf("%-60s ERROR fit-opt: %s | reset: %s\n", lab,
                substr(if (is.character(r1$p)) r1$p else "ok", 1, 90),
                substr(if (is.character(p2)) p2 else "ok", 1, 90)))
  } else {
    cat(sprintf("%-60s fit-opt vs in-sample %.3g ; reset vs fit-opt %.3g\n",
                lab, max(abs(unname(r1$p) - unname(r1$ins))),
                max(abs(unname(p2) - unname(r1$p)))))
  }
}
sumc <- c("contr.sum", "contr.poly")
trc <- c("contr.treatment", "contr.poly")
cmp("5353 fixed only y ~ f, sum then reset", y ~ f + (1 | g), sumc, trc)
cmp("5353 (1+f|g) sum then reset", y ~ f + (1 + f | g), sumc, trc)
cmp("5353 (1+o|g) ordered, poly->treatment for ordered",
    y ~ o + (1 + o | g), trc, c("contr.treatment", "contr.treatment"))
cmp("5353 (1+ch|g) character column", y ~ ch + (1 + ch | g), sumc, trc)
cmp("5353 (1+factor(cyl)|g) inline factor()", y ~ factor(cyl) +
      (1 + factor(cyl) | g), sumc, trc)
cmp("5353 (1+f:x|g) interaction", y ~ f * x + (1 + f:x | g), sumc, trc)
cmp("5353 (1+f|mm(g,g2))", y ~ x + (1 + f | mm(g, g2)), sumc, trc)
cmp("5353 (1+f||g)", y ~ f + (1 + f || g), sumc, trc)
cmp("5353 (1+f|gr(g, by=ch))", y ~ f + (1 + f | gr(g, by = ch)), sumc, trc)
cmp("5353 sigma ~ (1+f|g)", bf(y ~ f, sigma ~ (1 + f | g)), sumc, trc)
cmp("5353 helmert then reset", y ~ f + (1 + f | g),
    c("contr.helmert", "contr.poly"), trc)
# factor's own attribute at the fit, newdata without it
dca <- dc
contrasts(dca$f) <- contr.sum(3)
fa <- q(y ~ f + (1 + f | g), data = dca)
nd2 <- nd
nd2$f <- factor(as.character(nd2$f), levels = levels(dc$f))
pa <- tryCatch(fitted(fa, newdata = nd2)[, "Estimate"], error = function(e) conditionMessage(e))
cat(sprintf("%-60s %s\n", "5353 data attribute contr.sum, newdata plain",
            if (is.character(pa)) paste("ERROR", substr(pa, 1, 120)) else
              format(max(abs(unname(pa) - unname(fitted(fa)[c(1:6, 30), "Estimate"]))), digits = 3)))
# re_formula naming part of a term
withr::with_options(list(contrasts = sumc), {
  fp <- q(y ~ f + (1 + f | g) + (1 | g2), data = dc)
  r0 <- fitted(fp, newdata = nd, re_formula = ~ (1 + f | g))[, "Estimate"]
})
r1 <- tryCatch(fitted(fp, newdata = nd, re_formula = ~ (1 + f | g))[, "Estimate"],
               error = function(e) conditionMessage(e))
cat(sprintf("%-60s %s\n", "5353 re_formula partial after reset",
            if (is.character(r1)) paste("ERROR", substr(r1, 1, 120)) else
              format(max(abs(r1 - r0)), digits = 3)))
r1b <- tryCatch(fitted(fp, newdata = nd, re_formula = ~ (1 | g))[, "Estimate"],
               error = function(e) conditionMessage(e))
cat(sprintf("%-60s %s\n", "5353 re_formula (1|g) part of (1+f|g) after reset",
            if (is.character(r1b)) paste("ERROR", substr(r1b, 1, 120)) else "ran"))
cat("DONE\n")
