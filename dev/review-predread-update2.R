# Item 7, the rest of review-predread-update.R, whose print of a dropped-term prior
# table stopped it (a brms print defect, not a package one). Data seed 7601.
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
  if (inherits(obj, "error")) return(paste("ERR/REFUSED:", substr(conditionMessage(obj), 1, 220)))
  m <- as.matrix(obj)
  cn <- intersect(c("b_Intercept", "b_x", "b_z", "sigma"), colnames(m))
  paste(cn, round(colMeans(m[, cn, drop = FALSE]), 3), collapse = " ")
}
upb <- function(obj, ...) {
  tryCatch(quiet(update(obj, ..., refresh = 0, seed = 5, silent = 2)), error = function(e) e)
}
upd <- function(obj, ...) tryCatch(quiet(update(obj, ...)), error = function(e) e)
d0 <- samp(prior = sf("normal(10, 0.1)", class = "b", coef = "x") +
             sf("normal(0, 0.1)", class = "b", coef = "z"))
cat("=== formula change\n")
bj <- upb(bm, formula. = . ~ . - z)
cat("brms formula drop z:", ms(bj), "\n")
if (!inherits(bj, "error")) print(as.data.frame(bj$prior)[bj$prior$source == "user", c("prior", "class", "coef")])
dj <- upd(d0, formula. = . ~ . - z)
cat("draws formula drop z:", ms(dj), "\n")
dj2 <- upd(d0, formula. = . ~ . - z, prior = sf("normal(5, 0.1)", class = "sigma"))
cat("draws formula + prior:", ms(dj2), "\n")
cat("\n=== a new row naming a coef the model lacks\n")
bk <- upb(bm, prior = sb("normal(0,1)", class = "b", coef = "nosuch"))
cat("brms:", ms(bk), "\n")
dk <- upd(d0, prior = sf("normal(0,1)", class = "b", coef = "nosuch"))
cat("draws:", ms(dk), "\n")
dk0 <- tryCatch(samp(prior = sf("normal(0,1)", class = "b", coef = "nosuch")), error = function(e) e)
cat("frm_sample direct:", ms(dk0), "\n")
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
cat("\n=== sampling-route table: every row of it, as brms binds a default_prior table\n")
st <- tryCatch(frmtmb::default_prior(bf(y ~ x + z), family = gaussian(), data = d, route = "sample"),
               error = function(e) e)
cat("class:", if (inherits(st, "error")) conditionMessage(st) else paste(class(st), collapse = " "), " route:",
    if (!inherits(st, "error")) attr(st, "route"), "\n")
if (!inherits(st, "error")) {
  print(as.data.frame(st)[, c("prior", "class", "coef")])
  da <- upd(d0, prior = st)
  cat("frmtmb update(prior = whole table):", ms(da), "\n")
  tb <- brms::default_prior(brms::bf(y ~ x + z), data = d)
  bc <- upb(bm, prior = tb)
  cat("brms update(prior = whole table):", ms(bc), "\n")
  sub <- st[st$class == "b" & st$coef == "x", ]
  sub$prior <- "normal(-3, 0.1)"
  cat("class of a subset:", paste(class(sub), collapse = " "), " route:", attr(sub, "route"), "\n")
  dsub <- upd(d0, prior = sub)
  cat("frmtmb update(prior = one row of the table):", ms(dsub), "\n")
  tbs <- tb[tb$class == "b" & tb$coef == "x", ]
  tbs$prior <- "normal(-3, 0.1)"
  bsub <- upb(bm, prior = tbs)
  cat("brms update(prior = one row of the table):", ms(bsub), "\n")
  rt <- st
  rt$prior[rt$class == "b" & rt$coef == "x"] <- "normal(10, 0.1)"
  de <- tryCatch(samp(prior = rt), error = function(e) e)
  cat("first call with a table:", ms(de), "\n")
  if (!inherits(de, "error")) {
    df <- upd(de, prior = sf("normal(5, 0.1)", class = "sigma"))
    cat("table then sigma row:", ms(df), "\n")
    br <- tb
    br$prior[br$class == "b" & br$coef == "x"] <- "normal(10, 0.1)"
    bq <- quiet(brms::brm(brms::bf(y ~ x + z), data = d, chains = 1, iter = 300, refresh = 0, seed = 5, silent = 2, prior = br))
    bf1 <- upb(bq, prior = sb("normal(5, 0.1)", class = "sigma"))
    cat("brms table then sigma row:", ms(bf1), "\n")
    cat("brms table first:", ms(bq), "\n")
  }
}
cat("\n=== lb on a class\n")
dq <- tryCatch(samp(prior = sf("normal(10, 0.1)", class = "b", lb = 9.8)), error = function(e) e)
cat("class lb first:", ms(dq), "\n")
if (!inherits(dq, "error")) {
  dr <- upd(dq, prior = sf("normal(5, 0.1)", class = "sigma"))
  cat("then sigma:", ms(dr), "\n")
}
bq2 <- tryCatch(quiet(brms::brm(brms::bf(y ~ x + z), data = d, chains = 1, iter = 300, refresh = 0, seed = 5, silent = 2,
                                prior = sb("normal(10, 0.1)", class = "b", lb = 9.8))), error = function(e) e)
cat("brms class lb first:", ms(bq2), "\n")
if (!inherits(bq2, "error")) {
  br2 <- upb(bq2, prior = sb("normal(5, 0.1)", class = "sigma"))
  cat("brms then sigma:", ms(br2), "\n")
}
cat("\n=== rows differing in dpar under a sigma model\n")
ds0 <- quiet(frm_sample(bf(y ~ x, sigma ~ z), family = gaussian(), data = d, chains = 1, iter = 300, refresh = 0, seed = 5,
                        prior = sf("normal(0, 0.05)", class = "b", coef = "z", dpar = "sigma") +
                          sf("normal(10, 0.1)", class = "b", coef = "x")))
ds1 <- upd(ds0, prior = sf("normal(1, 0.05)", class = "b", coef = "z", dpar = "sigma"))
m0 <- as.matrix(ds1)
cat("after update: b_x", round(mean(m0[, "b_x"]), 3), " b_sigma_z", round(mean(m0[, "b_sigma_z"]), 3), "\n")
bs0 <- quiet(brms::brm(brms::bf(y ~ x, sigma ~ z), data = d, chains = 1, iter = 300, refresh = 0, seed = 5, silent = 2,
                       prior = sb("normal(0, 0.05)", class = "b", coef = "z", dpar = "sigma") +
                         sb("normal(10, 0.1)", class = "b", coef = "x")))
bs1 <- upb(bs0, prior = sb("normal(1, 0.05)", class = "b", coef = "z", dpar = "sigma"))
m1 <- as.matrix(bs1)
cat("brms after update: b_x", round(mean(m1[, "b_x"]), 3), " b_sigma_z", round(mean(m1[, "b_sigma_z"]), 3), "\n")
cat("\nfrmtmb attrs after merge:\n")
d1 <- upd(d0, prior = sf("normal(5, 0.1)", class = "sigma"))
sa <- d1$sample_args$prior
cat(" class", paste(class(sa), collapse = " "), " length", length(sa), " attrs",
    paste(names(attributes(sa)), collapse = " "), "\n")
print(tryCatch(frmtmb::prior_summary(d1), error = function(e) conditionMessage(e)))
