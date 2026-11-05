# Remaining single checks of attacks 1, 5, 8, 9.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages({library(rstan); library(brms)})
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
show <- function(label, expr) {
  r <- msg_of(eval.parent(substitute(expr)))
  cat("---", label, "\n")
  if (!is.null(r$err)) cat("ERR:", substr(r$err, 1, 300), "\n") else print(r$val)
  if (length(r$warn)) cat("WARN:", substr(paste(r$warn, collapse = " | "), 1, 240), "\n")
  invisible(r$val)
}
set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:12))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(12, arima.sim(list(ar = .5), 5, sd = .5)))
d$y <- rpois(n, exp(.5 + .4 * d$x + e))
d$bin <- rbinom(n, 1, plogis(.2 + .4 * d$x + e))
d$pos <- rgamma(n, 4, 4 / exp(.5 + .4 * d$x + e))
d$gs <- .4 * d$x + e + rnorm(n, 0, .5)
d$sub <- d$t != 3
d$ft <- factor(d$t, levels = rev(1:5))
d$ct <- letters[d$t]
d$dt <- d$t / 7
fp <- frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d)

cat("\n== coef() against brms on the same model\n")
pri <- c(set_prior("", class = "sderr"), set_prior("", class = "Intercept"))
form <- y ~ x + ar(t, g, cov = TRUE)
code <- brms::stancode(form, d, family = brms::brmsfamily("poisson"), prior = pri)
sdat <- brms::standata(form, d, family = brms::brmsfamily("poisson"), prior = pri)
m <- hel$brms_stan_model(code)
sf <- rstan::sampling(m, data = sdat, algorithm = "Fixed_param", chains = 1, iter = 20, warmup = 0, refresh = 0,
                      init = list(list(b = array(.4, 1), Intercept = .5, ar = array(.5, 1), zerr = rep(0, n), sderr = .5)))
scaf <- brm(form, d, family = brms::brmsfamily("poisson"), prior = pri, empty = TRUE)
scaf$fit <- sf
bf1 <- brms::rename_pars(scaf)
show("brms ranef", brms::ranef(bf1))
show("brms coef", brms::coef(bf1))
show("brms VarCorr", brms::VarCorr(bf1))
show("frm coef names", names(coef(fp)))
show("frm ngrps", ngrps(fp))

cat("\n== identity: factor time with reversed levels, decimal time, character time, duplicates\n")
mapf <- function(lab, form, dat, fam_b = brms::brmsfamily("poisson"), fam_f = poisson()) {
  br <- msg_of(brms::standata(form, dat, family = fam_b, prior = pri))
  fr <- msg_of(frm(form, data = dat, family = fam_f))
  cat(lab, ": brms", if (!is.null(br$err)) paste("ERR", substr(br$err, 1, 100)) else "ok", " frm", if (!is.null(fr$err)) paste("ERR", substr(fr$err, 1, 140)) else "ok", "\n")
  if (is.null(br$err)) {
    # the row order brms's data carries, by the response it was given
    cat("   brms Y first 8:", head(br$val$Y, 8), "\n")
  }
  if (is.null(fr$err) && is.null(br$err)) {
    a <- fi("frame_latres")(fr$val$frame)[[1]]$aux
    prev <- a$prev %||% a$lag_idx[[1]]
    ordr <- unlist(lapply(split(seq_len(nrow(dat)), a$group), function(i) {
      out <- integer(0)
      firsts <- i[prev[i] > nrow(dat) | prev[i] %in% (nrow(dat) + seq_len(a$G))]
      cur <- firsts[1]
      while (length(cur)) {
        out <- c(out, cur)
        cur <- i[prev[i] == cur]
      }
      out
    }))
    cat("   frm chain order Y  :", head(dat$y[ordr], 8), "\n")
    cat("   same order:", identical(as.numeric(br$val$Y), as.numeric(dat$y[ordr])), "\n")
  }
}
mapf("factor time, reversed levels", y ~ x + ar(ft, g, cov = TRUE), d)
mapf("decimal time", y ~ x + ar(dt, g, cov = TRUE), d)
mapf("character time", y ~ x + ar(ct, g, cov = TRUE), d)
du <- d
du$t[2] <- 1
mapf("duplicate time", y ~ x + ar(t, g, cov = TRUE), du)
set.seed(5)
ds <- d[sample(n), ]
mapf("unsorted, group as numeric", y ~ x + ar(t, gn, cov = TRUE), transform(ds, gn = as.integer(g)))
mapf("unsorted, gr = g, no time", y ~ x + ar(gr = g, cov = TRUE), ds)

cat("\n== mv with subset beside a latent; on_response osa; conditional_effects predict\n")
fm <- show("mv subset", {
  z <- frm(bf(y | subset(sub) ~ x + ar(t, g, cov = TRUE)) + bf(gs ~ x), data = d, family = poisson())
  logLik(z)
})
show("mv two latents", {
  z <- frm(bf(y ~ x + ar(t, g, cov = TRUE)) + bf(pos ~ x + cosy(t, g), family = Gamma(link = "log")), data = d, family = poisson())
  print(variables(z))
  logLik(z)
})
fg <- frm(pos ~ x + ar(t, g), family = Gamma(link = "log"), data = d)
show("osa on_response", summary(as.matrix(residuals(fg, type = "osa"))[, 1]))
fb <- frm(bin ~ x + ar(t, g, cov = TRUE), family = bernoulli(), data = d)
show("osa bernoulli", summary(as.matrix(residuals(fb, type = "osa"))[, 1]))
show("ce predict ML", head(conditional_effects(fp, "x", method = "posterior_predict")$x$estimate__, 2))
show("ce on_response", head(conditional_effects(fg, "x")$x$estimate__, 2))
show("importance plain", as.numeric(logLik(frm(y ~ x + (1 | g), family = poisson(), data = d, importance = 100))))
show("sd prior refused", frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d, prior = set_prior("normal(0, 1)", class = "sd")))
show("sderr prior moves", {
  f0 <- frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d, prior = set_prior("normal(0, 0.1)", class = "sderr"))
  c(prior = exp(f0$estimates$theta[1]), none = exp(fp$estimates$theta[1]))
})
show("sderr bound", {
  f0 <- frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d, prior = set_prior("", class = "sderr", lb = 0.9))
  exp(f0$estimates$theta[1])
})
show("simulate new missing t", simulate(fp, nsim = 1, newdata = d[1:3, c("x", "g")]))
cat("\n== frm_sample draws names, when a chain has more latents than rows\n")
suppressPackageStartupMessages(library(frmtmb.sample))
dr <- frm_sample(frmtmb::bf(y ~ x + ma(t, g, cov = TRUE)), family = poisson(), data = d, chains = 1, iter = 300, warmup = 200, seed = 11, refresh = 0)
v <- variables(dr)
cat("n rows", n, " n b[k]", length(grep("^b\\[", v)), "\n")
print(tail(v, 5))
