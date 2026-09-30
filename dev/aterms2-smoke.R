# Smoke test of cat(), rate(), subset(), index() and mi(idx =). Seed 7.
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
show <- function(label, expr) {
  cat("\n==========", label, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR:", conditionMessage(e)))
  print(r)
  invisible(r)
}
set.seed(7)
# cat
dc <- data.frame(s = sample(1:5, 40, TRUE), x = rnorm(40))
show("cat nthres", {
  fr <- frm(s | cat(6) ~ x, data = dc, family = cumulative(),
            dry_run = "frame")
  length(fr$par_template$tau_raw)
})
show("thres(5) nthres", {
  fr <- frm(s | thres(5) ~ x, data = dc, family = cumulative(),
            dry_run = "frame")
  length(fr$par_template$tau_raw)
})
show("cat + thres", frm(s | cat(6) + thres(5) ~ x, data = dc,
                        family = cumulative(), dry_run = "frame"))
show("cat gaussian", frm(s | cat(6) ~ x, data = dc, family = gaussian(),
                         dry_run = "frame"))

# rate
n <- 200
dr <- data.frame(x = rnorm(n), time = runif(n, 0.5, 3))
dr$y <- rpois(n, exp(0.3 + 0.5 * dr$x) * dr$time)
f1 <- frm(y | rate(time) ~ x, data = dr, family = poisson())
f2 <- frm(y ~ x + offset(log(time)), data = dr, family = poisson())
show("rate vs offset logLik", c(logLik(f1), logLik(f2),
                                 identical(as.numeric(logLik(f1)),
                                           as.numeric(logLik(f2)))))
show("rate vs offset coef", rbind(fixef(f1)[, 1], fixef(f2)[, 1]))
show("rate fitted vs offset fitted",
     max(abs(fitted(f1)[, 1] - fitted(f2)[, 1])))
show("rate fitted dpar mu / denom",
     range(fitted(f1)[, 1] / (fitted(f1, dpar = "mu")[, 1] * dr$time)))
show("rate predict newdata", fitted(f1, newdata = dr[1:3, ]))
show("rate predict newdata no time", fitted(f1, newdata = dr[1:3, "x", drop = FALSE]))
show("rate predict()", head(predict(f1)))
show("rate simulate", head(simulate(f1, nsim = 2, seed = 1)))
show("rate residuals", head(residuals(f1, type = "pearson")))
show("rate gaussian", frm(y | rate(time) ~ x, data = dr, family = gaussian()))
dr0 <- dr; dr0$time[1] <- 0
show("rate zero", frm(y | rate(time) ~ x, data = dr0, family = poisson()))
fnb <- frm(y | rate(time) ~ x, data = dr, family = negbinomial())
show("negbinomial rate", logLik(fnb))
fgeo <- frm(y | rate(time) ~ x, data = dr, family = geometric())
show("geometric rate", logLik(fgeo))
show("ce rate", head(conditional_effects(f1)[[1]]))

# subset, mv
set.seed(8)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(rep(letters[1:5], 12)),
                s1 = rep(c(TRUE, FALSE), 30), s2 = c(rep(TRUE, 45),
                                                     rep(FALSE, 15)))
d$y1 <- 1 + d$x + rnorm(5)[d$g] + rnorm(n)
d$y2 <- 2 - d$z + rnorm(n)
d$y2[!d$s2] <- NA
d$z[!d$s2] <- NA
bf1 <- bf(y1 | subset(s1) ~ x + (1 | g)) + bf(y2 | subset(s2) ~ z) +
  set_rescor(FALSE)
fs <- show("subset fit", frm(bf1, data = d, family = gaussian()))
show("subset frame", {
  fr <- frm(bf1, data = d, family = gaussian(), dry_run = "frame")
  c(n_obs = fr$n_obs, y1 = length(fr$y$y1), y2 = length(fr$y$y2),
    X_y2 = nrow(fr$linpreds[[2]]$X))
})
fa <- frm(y1 ~ x + (1 | g), data = d[d$s1, ], family = gaussian())
fb <- frm(y2 ~ z, data = d[d$s2, ], family = gaussian())
show("subset logLik vs separate", c(logLik(fs), logLik(fa) + logLik(fb)))
show("nobs", nobs(fs))
show("fitted all", fitted(fs))
show("fitted y1", dim(fitted(fs, resp = "y1")))
show("fitted y2", dim(fitted(fs, resp = "y2")))
show("predict y2", dim(predict(fs, resp = "y2")))
show("fitted newdata y1", fitted(fs, newdata = d[1:6, ], resp = "y1"))
show("summary", summary(fs))
show("ranef", ranef(fs))
show("fitted no resp", fitted(fs))
show("predict no resp", predict(fs))
show("frm_linpred no resp", frm_linpred(fs))
show("fitted newdata y2 no s1", dim(fitted(fs, newdata = d[1:50, setdiff(names(d), "s1")], resp = "y2")))
show("fitted newdata y1 no s1", fitted(fs, newdata = d[1:6, setdiff(names(d), "s1")], resp = "y1"))
show("predict newdata y1", dim(predict(fs, newdata = d[1:6, ], resp = "y1")))
show("ce", names(conditional_effects(fs)))
show("ce y2", head(conditional_effects(fs, effects = "z", resp = "y2")[[1]]))
show("VarCorr", VarCorr(fs))
show("confint", confint(fs))
show("residuals y1", dim(residuals(fs, resp = "y1")))
show("simulate", dim(simulate(fs)))
show("update", logLik(update(fs)))
show("vcov cluster", vcov_cluster(fs, cluster = ~ g))
show("rescor subset", frm(bf(y1 | subset(s1) ~ x) + bf(y2 | subset(s2) ~ z) +
                            set_rescor(TRUE), data = d, family = gaussian()))
d3 <- d; d3$s1[3] <- NA
show("NA subset", frm(bf1, data = d3, family = gaussian()))
show("univariate subset", {
  fu <- frm(y1 | subset(s1) ~ x, data = d, family = gaussian())
  fv <- frm(y1 ~ x, data = d[d$s1, ], family = gaussian())
  c(nobs(fu), logLik(fu), logLik(fv))
})
show("univariate subset newdata", fitted(frm(y1 | subset(s1) ~ x, data = d,
                                            family = gaussian()),
                                        newdata = d[1:6, ]))
# mi idx
set.seed(12)
dm <- data.frame(g1 = sample(seq(1, 59, 2), 60, TRUE), g2 = 1:60,
                 s = rep(c(TRUE, FALSE), 30))
dm$x <- rnorm(60)
dm$y <- 1 + 0.5 * dm$x[match(dm$g1, dm$g2)] + rnorm(60, sd = 0.5)
dm$x[c(3, 7, 8)] <- NA
bm <- bf(y ~ mi(x, idx = g1)) + bf(x | mi() + index(g2) + subset(s) ~ 1) +
  set_rescor(FALSE)
fm <- show("mi idx fit", frm(bm, data = dm, family = gaussian()))
show("mi idx frame", {
  fr <- frm(bm, data = dm, family = gaussian(), dry_run = "frame")
  list(idxl = fr$linpreds[[1]]$mi[[1]]$idxl, ny = length(fr$y$y),
       nx = length(fr$y$x), mi_map = fr$mi_map)
})
show("mi idx fitted y", dim(fitted(fm, resp = "y")))
show("mi idx fitted x", dim(fitted(fm, resp = "x")))
bm2 <- bf(y ~ mi(x, idx = g1)) + bf(x | mi() + subset(s) ~ 1) + set_rescor(FALSE)
show("no index", frm(bm2, data = dm, family = gaussian(), dry_run = "frame"))
bm3 <- bf(y ~ mi(x)) + bf(x | mi() + subset(s) + index(g2) ~ 1) + set_rescor(FALSE)
show("no idx, x subsetted", frm(bm3, data = dm, family = gaussian(), dry_run = "frame"))
bm4 <- bf(y | subset(s) ~ mi(x)) + bf(x | mi() ~ 1) + set_rescor(FALSE)
show("no idx, y subsetted", frm(bm4, data = dm, family = gaussian(), dry_run = "frame"))
dm5 <- dm; dm5$g1[2] <- 2
show("unmatched idx", frm(bm, data = dm5, family = gaussian(), dry_run = "frame"))
dm6 <- dm; dm6$g2[3] <- 1
show("duplicated index", frm(bm, data = dm6, family = gaussian(), dry_run = "frame"))
