# Review of lane refits: 3430 vcov_cluster() rows and update() overrides,
# 1567 dot formulas, 1259, 1285, 3427, 4375, 2939 names. Seeds inline.
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")

cat("\n==== 3430 vcov_cluster() with an outside cluster ====\n")
set.seed(1)
G <- 20
mkd <- function() {
  set.seed(1)
  dd <- data.frame(g = factor(rep(seq_len(G), each = 6)),
                   cl = factor(rep(seq_len(G / 2), each = 12)),
                   x = rnorm(G * 6))
  dd$y <- 1 + 0.5 * dd$x + rnorm(G * 6) + rnorm(G / 2)[dd$cl]
  dd
}
se_of <- function(v) round(sqrt(diag(v)), 7)
dd <- mkd()
fit <- frm(y ~ x, data = dd)
ref <- se_of(vcov_cluster(fit, dd$cl))
go("V0 explicit factor (reference)", ref)
go("V1 ~ cl, data present", se_of(vcov_cluster(fit, ~ cl)))
# NA rows dropped
dd <- mkd(); dd$x[c(3, 50, 51)] <- NA
fna <- frm(y ~ x, data = dd)
refna <- se_of(vcov_cluster(fna, dd$cl[-c(3, 50, 51)]))
go("V2 NA rows: ~ cl vs explicit", list(se_of(vcov_cluster(fna, ~ cl)), refna))
# reordered after the fit, row names kept
dd <- dd[sample(nrow(dd)), ]
go("V3 data reordered after fit, row names kept",
   list(se_of(vcov_cluster(fna, ~ cl)), refna))
# reordered, row names reset
rownames(dd) <- NULL
go("V4 data reordered, row names reset", se_of(vcov_cluster(fna, ~ cl)))
# subset in the call
dd <- mkd()
fsub <- frm(y ~ x, data = dd[dd$x > -1, ])
go("V5 data = dd[dd$x > -1, ]",
   list(se_of(vcov_cluster(fsub, ~ cl)),
        se_of(vcov_cluster(fsub, dd$cl[dd$x > -1]))))
# subset() term in the response
dd <- mkd(); dd$keep <- dd$x > -1
fsu <- go("V6 fit y | subset(keep) ~ x", frm(y | subset(keep) ~ x, data = dd))
if (!is.character(fsu)) {
  go("V6 ~ cl vs explicit",
     list(se_of(vcov_cluster(fsu, ~ cl)),
          se_of(vcov_cluster(fsu, dd$cl[dd$keep]))))
}
# tibble with NA rows
dd <- mkd(); dd$x[c(3, 50)] <- NA
tb <- tibble::as_tibble(dd)
ftb <- frm(y ~ x, data = tb)
go("V7 tibble NA rows: ~ cl vs explicit",
   list(se_of(vcov_cluster(ftb, ~ cl)), se_of(vcov_cluster(ftb, dd$cl[-c(3, 50)]))))
# character row names, NA rows
dd <- mkd(); rownames(dd) <- paste0("id", sample(1e4, nrow(dd)))
dd$x[c(7, 8)] <- NA
fch <- frm(y ~ x, data = dd)
go("V8 character row names: ~ cl vs explicit",
   list(se_of(vcov_cluster(fch, ~ cl)), se_of(vcov_cluster(fch, dd$cl[-c(7, 8)]))))
# the cluster variable changed after the fit (model variables untouched)
dd <- mkd(); fit <- frm(y ~ x, data = dd)
dd$cl <- factor(rep(seq_len(G / 4), each = 24))
go("V9 cl changed after fit (silent: not detectable)",
   list(se_of(vcov_cluster(fit, ~ cl)), ref))
# a global cl of length n beside the data's cl
dd <- mkd(); fit <- frm(y ~ x, data = dd)
cl <- factor(rep(1:2, length.out = nrow(dd)))
go("V10 global `cl` of length n beside data's cl",
   list(se_of(vcov_cluster(fit, ~ cl)), ref))
rm(cl)
# fit made in a function: data local, gone
mk <- function() { d_local <- mkd(); frm(y ~ x, data = d_local) }
fl <- mk()
go("V11 data local to a function, gone", se_of(vcov_cluster(fl, ~ cl)))
# fit made in a function, a global of the same name with other rows
d_local <- mkd()[1:60, ]
go("V12 data name resolves to another object", se_of(vcov_cluster(fl, ~ cl)))
rm(d_local)

cat("\n==== 3430 update() takes refit_args; explicit wins ====\n")
dd <- mkd(); dd$x[3] <- NA
mk_na <- function(d) { my_na <- stats::na.exclude; frm(y ~ x, data = d, na.action = my_na) }
fn <- mk_na(dd)
go("U1 update(fn) residual length (na.exclude kept)", length(residuals(update(fn))))
go("U2 update(fn, na.action = na.omit) residual length (explicit wins)",
   length(residuals(update(fn, na.action = stats::na.omit))))
mk_se <- function(d) { s <- TRUE; frm(y ~ x, data = d, se = s) }
fs <- mk_se(mkd())
go("U3 update(fs, se = FALSE)$se flag",
   { u <- update(fs, se = FALSE); c(is.null(u$sdr), getCall(u)$se) })
go("U4 update(fs) refit_args se", update(fs)$refit_args$se)
go("U5 update(fs, evaluate = FALSE)", update(fs, evaluate = FALSE))
go("U6 object size of refit_args (bytes)",
   c(fit = as.numeric(object.size(fs)),
     refit_args = as.numeric(object.size(fs$refit_args))))
go("U7 na.action closure env size in serialization",
   { mk_big <- function(d) { big <- rnorm(1e6); my_na <- function(object, ...) stats::na.omit(object, ...)
       frm(y ~ x, data = d, na.action = my_na) }
     fb <- mk_big(mkd())
     c(serialized_MB = length(serialize(fb, NULL)) / 1e6,
       refit_args_MB = length(serialize(fb$refit_args, NULL)) / 1e6) })

cat("\n==== 1567 dot formulas ====\n")
set.seed(3)
d3 <- data.frame(y = rnorm(40), x = rnorm(40), z = rnorm(40), w = rnorm(40),
                 o = abs(rnorm(40)))
d3n <- d3[, c("y", "x", "z", "o")]
cn <- function(f) rownames(fixef(f))
go("D1 y ~ . - w, with w", { f <- frm(y ~ . - w, data = d3); list(formula(f), cn(f)) })
go("D2 y ~ . - w, without w", cn(frm(y ~ . - w, data = d3n)))
go("D3 y ~ .^2 - w, with w", cn(frm(y ~ .^2 - w, data = d3)))
go("D4 y ~ .^2 - w, without w", cn(frm(y ~ .^2 - w, data = d3n)))
go("D5 update(y ~ .^2 - w, newdata without w)",
   cn(update(frm(y ~ .^2 - w, data = d3), newdata = d3n)))
go("D6 y ~ . - w - o + offset(o), with w",
   { f <- frm(y ~ . - w - o + offset(o), data = d3); list(formula(f), cn(f)) })
go("D7 same, update(newdata without w)",
   { f <- frm(y ~ . - w - o + offset(o), data = d3)
     u <- update(f, newdata = d3n)
     c(identical(fixef(u), fixef(f))) })
go("D8 y ~ . - q (q absent)", cn(frm(y ~ . - q, data = d3)))
go("D9 y ~ . - w + w:x", cn(frm(y ~ . - w + w:x, data = d3)))
go("D10 log(y + 5) ~ . - w", cn(frm(log(y + 5) ~ . - w, data = d3)))
go("D11 y ~ -w + .", cn(frm(y ~ -w + ., data = d3)))
go("D12 y ~ . - w - x - z - o (all removed)", { f <- frm(y ~ . - w - x - z - o, data = d3); list(formula(f), cn(f)) })
go("D13 y ~ . - w in sigma", { f <- frm(bf(y ~ x, sigma ~ . - w - x - y), data = d3); list(cn(f)) })
go("D14 y ~ 0 + . - w", cn(frm(y ~ 0 + . - w, data = d3)))
go("D15 y ~ . - w, w a factor, update without w",
   { d4 <- d3; d4$w <- factor(d4$w > 0); f <- frm(y ~ . - w, data = d4)
     cn(update(f, newdata = d4[, -4])) })

cat("\n==== 1259 character ordinal ====\n")
set.seed(2501)
d1 <- data.frame(x = rnorm(60))
d1$yc <- sample(c("none", "mild", "moderate", "severe"), 60, TRUE)
d1$yd <- sample(c("1", "2", "3", "4"), 60, TRUE)
go("O1 labels", frm(yc ~ x, data = d1, family = cumulative()))
go("O2 digits", round(fixef(frm(yd ~ x, data = d1, family = cumulative()))[, 1], 6))
go("O3 digits vs integers", {
  a <- fixef(frm(yd ~ x, data = d1, family = cumulative()))[, 1]
  b <- fixef(frm(as.integer(yd) ~ x, data = transform(d1, yd = as.integer(yd)),
                 family = cumulative()))[, 1]
  max(abs(a - b)) })
go("O4 digits with NA and 'n/a'", {
  d1b <- d1; d1b$yd[1:3] <- c(NA, "n/a", "")
  frm(yd ~ x, data = d1b, family = cumulative()) })
go("O5 decimal text", { d1c <- d1; d1c$yd[1] <- "2.5"; frm(yd ~ x, data = d1c, family = cumulative()) })
go("O6 hurdle_cumulative labels? (code0 = 0 wording)", {
  frm(yc ~ x, data = d1, family = sratio()) })
go("O7 categorical labels unaffected", class(frm(yc ~ x, data = d1, family = categorical())))

cat("\n==== 1285 fitted mv ordinal one row ====\n")
set.seed(77)
n <- 300
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$yo <- sample(1:3, n, TRUE); d2$yo2 <- sample(1:3, n, TRUE)
mv2 <- suppressWarnings(frm(bf(yo ~ x) + bf(yo2 ~ z), family = sratio(), data = d2))
nd1 <- data.frame(x = 0.3, z = -0.2)
nd2 <- data.frame(x = c(0.3, 1), z = c(-0.2, 0))
go("M1 fitted 1 row == row 1 of 2 rows; dims",
   { a <- fitted(mv2, newdata = nd1); b <- fitted(mv2, newdata = nd2)
     list(dim(a), dimnames(a)[[2]], max(abs(a[1, , ] - b[1, , ]))) })
go("M2 fitted(summary = FALSE) 1 row", dim(fitted(mv2, newdata = nd1, summary = FALSE)))
go("M3 predict 1 row", dim(predict(mv2, newdata = nd1)))

cat("\n==== 3427 getCall() after a delta update ====\n")
set.seed(7)
dk <- data.frame(x = rnorm(40), x3 = rnorm(40), g = factor(rep(1:8, 5)))
dk$z <- dk$x + rnorm(8)[dk$g] + rnorm(40)
fk <- frm(z ~ x + x3 + (1 | g), data = dk, family = student())
go("C1 family arg student, delta", { u <- update(fk, ~ . - x3); print(getCall(u))
  e <- eval(getCall(u)); identical(fixef(e), fixef(u)) })
fk2 <- frm(bf(z ~ x + x3 + (1 | g)) + student(), data = dk)
go("C2 family attached in bf, delta", { u <- update(fk2, ~ . - x3); print(getCall(u))
  e <- eval(getCall(u)); c(identical(fixef(e), fixef(u)), family(u)$family) })
go("C3 delta with a new family", { u <- update(fk2, ~ . - x3, family = gaussian()); print(getCall(u))
  e <- eval(getCall(u)); c(identical(fixef(e), fixef(u)), family(u)$family) })
go("C4 delta on bf with sigma ~ and autocor", {
  fs <- frm(bf(z ~ x + x3, sigma ~ x, center = FALSE), data = dk)
  u <- update(fs, ~ . - x3); print(getCall(u))
  e <- eval(getCall(u)); identical(fixef(e), fixef(u)) })

cat("\n==== 4375 saved ar(cov = TRUE) ====\n")
set.seed(4)
dh <- data.frame(t = rep(1:12, 5), h = factor(rep(1:5, each = 12)), x = rnorm(60))
dh$y <- 1 + 0.5 * dh$x + rnorm(60)
fh <- frm(y ~ x + ar(time = t, gr = h, cov = TRUE), data = dh)
old <- fh; old$frame$autocor[[1L]]$sigma_is <- NULL
go("A1 confint(old, parm = 'x')", confint(old, parm = "x"))
go("A2 confint(old, parm = 'x', method = 'profile')", confint(old, parm = "x", method = "profile"))
go("A3 confint(fh, parm = 'x', method = 'profile')", confint(fh, parm = "x", method = "profile"))
go("A4 confint(old, parm = 'x', method = 'uniroot')", confint(old, parm = "x", method = "uniroot"))
go("A5 confint(old, parm = 'sigma')", confint(old, parm = "sigma"))
go("A6 confint(old, parm = 2)", confint(old, parm = 2))
go("A7 parnames(old)", parnames(old))
go("A8 confint(old) with every index but sigma",
   { pn <- rownames(confint(fh)); confint(old, parm = setdiff(pn, grep("sigma", pn, value = TRUE))) })
go("A9 confint(old, parm = 'x', vcov = vcov_cluster(old, ~ h))",
   confint(old, parm = "x", vcov = vcov_cluster(old, ~ h)))
fh2 <- frm(bf(y ~ x + ar(time = t, gr = h, cov = TRUE), sigma ~ x), data = dh)
old2 <- fh2; old2$frame$autocor[[1L]]$sigma_is <- NULL
go("A10 sigma ~ x: confint(old2, parm = 'x')", confint(old2, parm = "x"))
go("A11 sigma ~ x: confint(old2, parm = 'sigma_x')", confint(old2, parm = "sigma_x"))
go("A12 parnames(old2)", parnames(old2))
