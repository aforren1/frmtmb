# Item 7: update(<draws>, prior = ) against brms's update.R 144-155. Data seed 7601.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7601)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, 20)))
d$y <- 10 + 2 * d$x - d$z + rnorm(n)
sb <- brms::set_prior
sf <- frmtmb::set_prior
samp <- function(...) {
  quiet(frm_sample(bf(y ~ x + z), family = gaussian(), data = d, chains = 1,
                   iter = 300, refresh = 0, seed = 5, ...))
}
bm <- quiet(brms::brm(brms::bf(y ~ x + z), data = d, chains = 1, iter = 300, refresh = 0,
                      seed = 5, silent = 2,
                      prior = sb("normal(10, 0.1)", class = "b", coef = "x") +
                        sb("normal(0, 0.1)", class = "b", coef = "z")))
ms <- function(obj) {
  if (inherits(obj, "error")) return(paste("ERR/REFUSED:", substr(conditionMessage(obj), 1, 200)))
  m <- as.matrix(obj)
  cn <- intersect(c("b_Intercept", "b_x", "b_z", "sigma"), colnames(m))
  paste(cn, round(colMeans(m[, cn, drop = FALSE]), 3), collapse = " ")
}
pr <- function(b) {
  p <- b$prior
  p[p$source == "user", c("prior", "class", "coef", "source")]
}
upb <- function(obj, ...) {
  tryCatch(quiet(update(obj, ..., refresh = 0, seed = 5, silent = 2)), error = function(e) e)
}
upd <- function(obj, ...) tryCatch(quiet(update(obj, ...)), error = function(e) e)
cat("=== base: b_x ~ N(10, .1), b_z ~ N(0, .1)\n")
d0 <- samp(prior = sf("normal(10, 0.1)", class = "b", coef = "x") +
             sf("normal(0, 0.1)", class = "b", coef = "z"))
cat("brms  ", ms(bm), "\nfrmtmb", ms(d0), "\n")
cat("\n=== new row in another slot: sigma\n")
b1 <- upb(bm, prior = sb("normal(5, 0.1)", class = "sigma"))
d1 <- upd(d0, prior = sf("normal(5, 0.1)", class = "sigma"))
cat("brms  ", ms(b1), "\nfrmtmb", ms(d1), "\n")
print(pr(b1))
cat("\n=== replace a row in its slot: b_x ~ N(-3, .1)\n")
b2 <- upb(bm, prior = sb("normal(-3, 0.1)", class = "b", coef = "x"))
d2 <- upd(d0, prior = sf("normal(-3, 0.1)", class = "b", coef = "x"))
cat("brms  ", ms(b2), "\nfrmtmb", ms(d2), "\n")
print(pr(b2))
cat("\n=== duplicate rows in new prior\n")
b3 <- upb(bm, prior = sb("normal(-3, 0.1)", class = "b", coef = "x") +
            sb("normal(-6, 0.1)", class = "b", coef = "x"))
d3 <- upd(d0, prior = sf("normal(-3, 0.1)", class = "b", coef = "x") +
            sf("normal(-6, 0.1)", class = "b", coef = "x"))
cat("brms  ", ms(b3), "\nfrmtmb", ms(d3), "\n")
cat("\n=== class-wide new row beside old coef rows: class b ~ N(50,.1)\n")
b4 <- upb(bm, prior = sb("normal(50, 0.1)", class = "b"))
d4 <- upd(d0, prior = sf("normal(50, 0.1)", class = "b"))
cat("brms  ", ms(b4), "\nfrmtmb", ms(d4), "\n")
print(pr(b4))
cat("\n=== remove: prior '' on b_x\n")
b5 <- upb(bm, prior = sb("", class = "b", coef = "x"))
d5 <- upd(d0, prior = sf("", class = "b", coef = "x"))
cat("brms  ", ms(b5), "\nfrmtmb", ms(d5), "\n")
if (!inherits(b5, "error")) print(pr(b5))
cat("\n=== prior = NULL\n")
b6 <- upb(bm, prior = NULL)
d6 <- upd(d0, prior = NULL)
cat("brms  ", ms(b6), "\nfrmtmb", ms(d6), "\n")
cat("\n=== empty_prior()\n")
b7 <- upb(bm, prior = brms::empty_prior())
d7 <- upd(d0, prior = frmtmb::empty_prior())
cat("brms  ", ms(b7), "\nfrmtmb", ms(d7), "\n")
cat("\n=== sequential updates: sigma, then b_x\n")
b8 <- upb(b1, prior = sb("normal(-3, 0.1)", class = "b", coef = "x"))
d8 <- upd(d1, prior = sf("normal(-3, 0.1)", class = "b", coef = "x"))
cat("brms  ", ms(b8), "\nfrmtmb", ms(d8), "\n")
cat("\n=== sample_prior kept; prior draws under the new prior\n")
bsp <- quiet(brms::brm(brms::bf(y ~ x + z), data = d, chains = 1, iter = 300, refresh = 0,
                       seed = 5, silent = 2, sample_prior = "yes",
                       prior = sb("normal(10, 0.1)", class = "b", coef = "x") +
                         sb("normal(0, 1)", class = "b", coef = "z")))
dsp <- samp(sample_prior = "yes", prior = sf("normal(10, 0.1)", class = "b", coef = "x") +
              sf("normal(0, 1)", class = "b", coef = "z"))
b9 <- upb(bsp, prior = sb("normal(-3, 0.1)", class = "b", coef = "x"))
d9 <- upd(dsp, prior = sf("normal(-3, 0.1)", class = "b", coef = "x"))
cat("brms prior_ columns:", paste(grep("^prior_", colnames(as.matrix(b9)), value = TRUE), collapse = " "), "\n")
pd <- tryCatch(frmtmb.sample::prior_draws(d9), error = function(e) e)
cat("draws prior_draws:", if (inherits(pd, "error")) conditionMessage(pd) else paste(colnames(pd), collapse = " "), "\n")
cat("prior_b_x mean: brms ", round(mean(as.matrix(b9)[, "prior_b_x"]), 3), " draws ",
    if (inherits(pd, "error")) "ERR" else round(colMeans(pd[, grep("b_x", colnames(pd)), drop = FALSE]), 3), "\n")
cat("prior_b_z mean: brms ", round(mean(as.matrix(b9)[, "prior_b_z"]), 3), " sd ", round(sd(as.matrix(b9)[, "prior_b_z"]), 3), " draws sd ",
    if (inherits(pd, "error")) "ERR" else round(apply(pd[, grep("b_z", colnames(pd)), drop = FALSE], 2, sd), 3), "\n")
d10 <- upd(d0, sample_prior = "yes")
b10 <- upb(bm, sample_prior = "yes")
cat("sample_prior yes only: brms ", ms(b10), " draws ", ms(d10), "\n")
cat("\n=== route marks: a fit-route table\n")
ft <- frmtmb::default_prior(bf(y ~ x + z), family = gaussian(), data = d)
cat("class of core default_prior:", paste(class(ft), collapse = " "), " route attr:",
    paste(names(attributes(ft)), collapse = " "), "\n")
da <- tryCatch(samp(prior = ft), error = function(e) e)
cat("frm_sample(prior = core table):", ms(da), "\n")
db <- upd(d0, prior = ft)
cat("update(<draws>, prior = core table):", ms(db), "\n")
st <- tryCatch(frmtmb.sample::default_prior(bf(y ~ x + z), family = gaussian(), data = d, sample = TRUE),
               error = function(e) e)
if (inherits(st, "error")) {
  st <- tryCatch(frmtmb.sample::sample_default_prior(bf(y ~ x + z), family = gaussian(), data = d),
                 error = function(e) e)
}
cat("sampling-route table:", if (inherits(st, "error")) conditionMessage(st) else paste(class(st), collapse = " "), "\n")
tb <- brms::default_prior(brms::bf(y ~ x + z), data = d)
bc <- upb(bm, prior = tb)
cat("brms update(prior = default_prior table):", ms(bc), "\n")
cat("\n=== the fit's own prior (frm(prior =)), then frm_sample(prior =), then update(prior =)\n")
fo <- quiet(frm(bf(y ~ x + z), family = gaussian(), data = d,
                prior = sf("normal(10, 0.1)", class = "b", coef = "x")))
dg <- quiet(frm_sample(fo, chains = 1, iter = 300, refresh = 0, seed = 5,
                       prior = sf("normal(0, 0.1)", class = "b", coef = "z")))
cat("own + sample:", ms(dg), "\n")
dh <- upd(dg, prior = sf("normal(5, 0.1)", class = "sigma"))
cat("own + sample + update sigma:", ms(dh), "\n")
di <- upd(dg, prior = sf("normal(-3, 0.1)", class = "b", coef = "x"))
cat("update on an own slot:", ms(di), "\n")
cat("\n=== formula change\n")
bj <- upb(bm, formula. = . ~ . - z)
cat("brms formula drop z:", ms(bj), "\n")
if (!inherits(bj, "error")) print(pr(bj))
dj <- upd(d0, formula. = . ~ . - z)
cat("draws formula drop z:", ms(dj), "\n")
cat("\n=== a new row naming a coef the model lacks\n")
bk <- upb(bm, prior = sb("normal(0,1)", class = "b", coef = "nosuch"))
cat("brms:", ms(bk), "\n")
dk <- upd(d0, prior = sf("normal(0,1)", class = "b", coef = "nosuch"))
cat("draws:", ms(dk), "\n")
cat("\n=== a first call with no prior at all\n")
dn <- samp()
dm <- upd(dn, prior = sf("normal(-3, 0.1)", class = "b", coef = "x"))
cat("no prior then b_x:", ms(dm), "\n")
cat("\n=== update(<draws>) with no prior argument\n")
do <- upd(d0, iter = 300)
cat("no prior arg:", ms(do), "\n")
cat("\n=== a non-brmsprior: brms refuses\n")
bl <- upb(bm, prior = "normal(0,1)")
cat("brms:", ms(bl), "\n")
dl <- upd(d0, prior = "normal(0,1)")
cat("draws:", ms(dl), "\n")
cat("\n=== lb on a coef\n")
dq <- tryCatch(samp(prior = sf("normal(10, 0.1)", class = "b", lb = 9.8)), error = function(e) e)
cat("class lb first:", ms(dq), "\n")
if (!inherits(dq, "error")) {
  dr <- upd(dq, prior = sf("normal(5, 0.1)", class = "sigma"))
  cat("then sigma:", ms(dr), "\n")
}
cat("\nfrmtmb attrs after merge:\n")
sa <- d1$sample_args$prior
cat(" class", paste(class(sa), collapse = " "), " length", length(sa), " attrs",
    paste(names(attributes(sa)), collapse = " "), "\n")
print(tryCatch(frmtmb::prior_summary(d1), error = function(e) conditionMessage(e)))
