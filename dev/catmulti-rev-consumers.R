# Reviewer: every consumer of a non-first refcat and of refcat = NA, on
# the lane library. Log dev/catmulti-rev-out/consumers.txt.
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
source("dev/catmulti-rev-data.R")
sec <- function(t) cat("\n========", t, "========\n")
shape <- function(nm, x) {
  if (is.null(dim(x))) {
    cat(nm, ": length", length(x), "class", class(x)[1], "\n")
    if (!is.null(names(x))) cat("  names:", head(names(x), 8), "\n")
    return(invisible())
  }
  cat(nm, ": dim", paste(dim(x), collapse = " x "), "\n")
  dn <- dimnames(x)
  for (k in seq_along(dn)) {
    cat("  dimnames[[", k, "]]: ",
        if (is.null(dn[[k]])) "NULL" else
          paste(utils::head(dn[[k]], 8), collapse = ", "), "\n", sep = "")
  }
}
TRY <- function(nm, expr) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "err")),
    warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") },
    message = function(cm) invokeRestart("muffleMessage"))
  if (inherits(r, "err")) cat(nm, ": ERROR:", r, "\n") else if (is.list(r) &&
      !is.data.frame(r)) { cat(nm, ": list", names(r), "\n") } else shape(nm, r)
  if (length(w)) cat("  WARNINGS:", substr(unique(w), 1, 200), sep = "\n    ")
  invisible(if (inherits(r, "err")) NULL else r)
}
rel <- function(a, b) max(abs(a - b)) / max(abs(a), abs(b))

d3 <- rev_cat_data(21, 3)
d4 <- rev_cat_data(22, 4)
newg <- data.frame(x = c(-1, 0, 1), g = c("g1", "g2", "new"))

for (K in c(3, 4)) {
  dd <- if (K == 3) d3 else d4
  sec(paste0("K = ", K, " (1 | p | g), refcat = 'b' vs default"))
  fb <- TRY("fit refb", frm(bf(y ~ x + (1 | p | g)),
                           family = categorical(refcat = "b"), data = dd))
  fa <- TRY("fit refa", frm(bf(y ~ x + (1 | p | g)),
                           family = categorical(), data = dd))
  cat("logLik refb", format(as.numeric(logLik(fb)), digits = 12),
      "refa", format(as.numeric(logLik(fa)), digits = 12), "\n")
  cat("variables refb:\n"); print(grep("^(b_|sd_|cor_)", variables(fb), value = TRUE))
  print(round(summary(fb)$fixed, 4))
  print(VarCorr(fb))
  rb <- TRY("ranef refb", ranef(fb)$g)
  cb <- TRY("coef refb", coef(fb)$g)
  if (!is.null(cb)) print(round(cb[1:2, , , drop = FALSE], 4))
  Pb <- TRY("fitted refb", fitted(fb))
  Pa <- fitted(fa)
  cat("max |P refb - P refa| (with RE):", max(abs(Pb[, "Estimate", ] - Pa[, "Estimate", ])), "\n")
  cat("max |SE refb - SE refa|/SE:", max(abs(Pb[, "Est.Error", ] - Pa[, "Est.Error", ]) /
                                       Pa[, "Est.Error", ]), "\n")
  P0b <- fitted(fb, re_formula = NA); P0a <- fitted(fa, re_formula = NA)
  cat("max |P refb - P refa| (re_formula NA):",
      max(abs(P0b[, "Estimate", ] - P0a[, "Estimate", ])), "\n")
  TRY("fitted linear refb", fitted(fb, scale = "linear"))
  TRY("predict refb", predict(fb))
  pr <- TRY("predict refb head", predict(fb)[1:3, , drop = FALSE])
  if (!is.null(pr)) print(pr)
  set.seed(5)
  s <- TRY("simulate refb", simulate(fb, nsim = 200, seed = 5))
  if (!is.null(s)) {
    tab <- table(factor(unlist(s), levels = levels(dd$y))) / (200 * nrow(dd))
    cat("simulate shares:", round(tab, 4), "\n")
    cat("fitted mean P  :", round(colMeans(Pb[, "Estimate", ]), 4), "\n")
    cat("observed shares:", round(prop.table(table(dd$y)), 4), "\n")
  }
  ce <- TRY("conditional_effects", conditional_effects(fb, "x", categorical = TRUE)[[1]])
  if (!is.null(ce)) { cat("ce cols:", names(ce), "\n"); cat("ce cats:", levels(ce$cats__), "\n") }
  nf <- TRY("fitted newdata new level", fitted(fb, newdata = newg, allow_new_levels = TRUE))
  nf0 <- TRY("fitted newdata re NA", fitted(fb, newdata = newg, re_formula = NA))
  if (!is.null(nf0)) print(round(nf0[, "Estimate", ], 5))
  TRY("fitted newdata new level refused", fitted(fb, newdata = newg))
  TRY("predict newdata new level", predict(fb, newdata = newg, allow_new_levels = TRUE))
  TRY("frm_linpred dpar mua newdata", frm_linpred(fb, newdata = newg, dpar = "mua",
                                                  re_formula = NA))
  TRY("frm_linpred dpar mub (reference) refused", frm_linpred(fb, dpar = "mub"))
  if (requireNamespace("emmeans", quietly = TRUE)) {
    TRY("emmeans dpar mua", summary(emmeans::emmeans(fb, ~ x, dpar = "mua", at = list(x = 0))))
    em <- TRY("emmeans default", summary(emmeans::emmeans(fb, ~ x, at = list(x = 0))))
    if (!is.null(em)) print(em)
    em2 <- TRY("emmeans epred", summary(emmeans::emmeans(fb, ~ x, epred = TRUE, at = list(x = 0))))
    if (!is.null(em2)) print(em2)
  }
  fu <- TRY("update refcat c", update(fb, family = categorical(refcat = "c")))
  if (!is.null(fu)) {
    cat("update dpars:", names(fixef_by_dpar(fu)), " logLik", format(as.numeric(logLik(fu)), digits = 12), "\n")
  }
  fu2 <- TRY("update data", update(fb, data = dd[1:300, ]))
  if (!is.null(fu2)) cat("update(data) keeps refcat:", family(fu2)$refcat,
                         "dpars", names(fixef_by_dpar(fu2)), "\n")
  fu3 <- TRY("update formula", update(fb, . ~ . - x))
  if (!is.null(fu3)) cat("update(formula) refcat:", family(fu3)$refcat,
                         "dpars", names(fixef_by_dpar(fu3)), "\n")
  TRY("default_prior refb", print(default_prior(bf(y ~ x + (1 | p | g)), data = dd,
                                                family = categorical(refcat = "b"))))
  TRY("drop1 refb", print(drop1(fb, test = "Chisq")))
}

sec("refcat = NA, (1 | g)")
fn <- TRY("fit refNA RE", frm(bf(y ~ x + (1 | p | g)), family = categorical(refcat = NA), data = d3))
if (!is.null(fn)) {
  TRY("fitted refNA", fitted(fn)); TRY("fitted linear refNA", fitted(fn, scale = "linear"))
  TRY("ranef refNA", ranef(fn)$g)
  fa <- frm(bf(y ~ x + (1 | p | g)), family = categorical(), data = d3)
  cat("logLik NA", as.numeric(logLik(fn)), "ref a", as.numeric(logLik(fa)), "\n")
}

sec("own formula, reference dpar, priors under refcat = 'b'")
TRY("bf(y ~ x, muc ~ 1)", fixef_by_dpar(frm(bf(y ~ x, muc ~ 1), family = categorical(refcat = "b"), data = d3)))
TRY("bf(y ~ x, mub ~ 1) names the reference", frm(bf(y ~ x, mub ~ 1), family = categorical(refcat = "b"), data = d3))
TRY("prior on reference dpar", frm(bf(y ~ x), family = categorical(refcat = "b"), data = d3,
                                   prior = set_prior("normal(0, 1)", class = "b", dpar = "mub")))
fp <- TRY("prior on mua", frm(bf(y ~ x), family = categorical(refcat = "b"), data = d3,
                              prior = set_prior("normal(0, 0.1)", class = "b", dpar = "mua")))
if (!is.null(fp)) print(round(fixef(fp), 4))

sec("mixture / hurdle with categorical")
TRY("mixture(categorical, categorical)", mixture(categorical(), categorical()))
TRY("mixture(multinomial, multinomial)", mixture(multinomial(), multinomial()))

sec("make_stan_names edge cases")
for (lv in list(c("a.b", "ab", "c"), c("a b", "1x", "c"), c("a_b", "c.d", "e"),
                c("x", "x.1", "x1"))) {
  set.seed(3)
  d <- data.frame(x = rnorm(60))
  d$y <- factor(sample(lv, 60, TRUE), levels = lv)
  r <- TRY(paste("levels", paste(lv, collapse = " | ")),
           names(fixef_by_dpar(frm(bf(y ~ x), family = categorical(), data = d))))
  if (!is.null(r)) cat("  ->", r, "\n")
  r <- TRY("  refcat NA", names(fixef_by_dpar(suppressWarnings(
    frm(bf(y ~ x), family = categorical(refcat = NA), data = d)))))
  if (!is.null(r)) cat("  ->", r, "\n")
  r <- TRY("  default_prior", unique(default_prior(bf(y ~ x), data = d, family = categorical())$dpar))
  if (!is.null(r)) cat("  ->", r, "\n")
}
set.seed(3)
d <- data.frame(x = rnorm(60)); d$y <- factor(sample(c("", "a", "b"), 60, TRUE))
r <- TRY("level ''", names(fixef_by_dpar(frm(bf(y ~ x), family = categorical(refcat = "a"), data = d))))
if (!is.null(r)) cat("  ->", r, "\n")

sec("bare multinomial()")
dm <- rev_mult_data(31)
for (f in c("cbind(y1, y2, y3) | trials(N) ~ x",
            "cbind(y1, y2 + 0, y3) | trials(N) ~ x",
            "cbind(a = y1, `b c` = y2, y3) | trials(N) ~ x",
            "Yc | trials(N) ~ x")) {
  r <- TRY(f, names(fixef_by_dpar(frm(bf(as.formula(f)), family = multinomial(), data = dm))))
  if (!is.null(r)) cat("  ->", r, "\n")
}
dn <- data.frame(x = dm$x, N = dm$N, check.names = FALSE)
dn[["y 1"]] <- dm$y1; dn[["y-2"]] <- dm$y2; dn[["3y"]] <- dm$y3
r <- TRY("non-syntactic names", names(fixef_by_dpar(frm(bf(cbind(`y 1`, `y-2`, `3y`) | trials(N) ~ x),
                                                           family = multinomial(), data = dn))))
if (!is.null(r)) cat("  ->", r, "\n")
dm2 <- dm; dm2$Yu <- unname(dm$Yc); colnames(dm2$Yu) <- NULL
r <- TRY("unnamed matrix column", names(fixef_by_dpar(frm(bf(Yu | trials(N) ~ x), family = multinomial(), data = dm2))))
if (!is.null(r)) cat("  ->", r, "\n")
r <- TRY("multinomial(refcat = 2) unnamed", names(fixef_by_dpar(frm(bf(Yu | trials(N) ~ x),
                                                                    family = multinomial(refcat = 2), data = dm2))))
if (!is.null(r)) cat("  ->", r, "\n")
TRY("multinomial(K = 4) on 3 columns", frm(bf(Yc | trials(N) ~ x), family = multinomial(K = 4), data = dm))
r <- TRY("multinomial(K = 3) on named columns", names(fixef_by_dpar(frm(bf(Yc | trials(N) ~ x), family = multinomial(K = 3), data = dm))))
if (!is.null(r)) cat("  ->", r, "\n")
TRY("multinomial(3)", multinomial(3))
TRY("multinomial('logit')", multinomial("logit")$family)
TRY("multinomial(logit)", multinomial(logit)$family)
TRY("multinomial('probit')", multinomial("probit"))
TRY("brmsfamily('multinomial', refcat = 'y2')", brmsfamily("multinomial", refcat = "y2")$refcat)
TRY("brms::multinomial(refcat = 'y3') fit", names(fixef_by_dpar(frm(bf(Yc | trials(N) ~ x),
                                                                    family = brms::multinomial(refcat = "y3"), data = dm))))
TRY("no trials", frm(bf(Yc ~ x), family = multinomial(), data = dm))
TRY("one column", frm(bf(cbind(y1) | trials(N) ~ x), family = multinomial(), data = dm))
TRY("rows not summing to trials", frm(bf(Yc | trials(N + 1) ~ x), family = multinomial(), data = dm))

sec("multinomial fitted, varying trials")
fm <- frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dm)
print(round(fixef(fm), 5))
ft <- fitted(fm)
shape("fitted", ft)
print(round(ft[1:6, , 1], 4)); print(round(ft[1:6, , 3], 4))
cat("trials of rows 1..6:", dm$N[1:6], "\n")
i1 <- which(dm$N == 1)[1:3]
print(round(ft[i1, , 3], 5))
# independent delta method: E = N softmax(X B), SE from vcov
V <- vcov(fm)
cf <- fixef(fm)[, "Estimate"]
cat("coef names:", names(cf), "\n")
cat("vcov names:", rownames(V), "\n")
X <- cbind(1, dm$x)
idx2 <- c(grep("muy2_Intercept", rownames(V)), grep("muy2_x", rownames(V)))
idx3 <- c(grep("muy3_Intercept", rownames(V)), grep("muy3_x", rownames(V)))
B2 <- cf[idx2]; B3 <- cf[idx3]
e2 <- drop(X %*% B2); e3 <- drop(X %*% B3)
den <- 1 + exp(e2) + exp(e3)
P <- cbind(1, exp(e2), exp(e3)) / den
E <- P * dm$N
se <- matrix(NA, nrow(dm), 3)
Vs <- V[c(idx2, idx3), c(idx2, idx3)]
for (i in seq_len(nrow(dm))) {
  p <- P[i, ]
  # dP_k / d eta_j = p_k (1{k=j} - p_j), j = 2, 3
  dPde <- cbind(p * ((1:3 == 2) - p[2]), p * ((1:3 == 3) - p[3]))
  J <- dm$N[i] * cbind(dPde[, 1] %o% X[i, ], dPde[, 2] %o% X[i, ])
  J <- matrix(J, 3)
  se[i, ] <- sqrt(diag(J %*% Vs %*% t(J)))
}
cat("estimate rel diff vs independent:", rel(E, ft[, "Estimate", ]), "\n")
cat("SE rel diff vs independent (max over cells):",
    max(abs(se - ft[, "Est.Error", ]) / se), "\n")
lo <- ft[, "Q2.5", ]; hi <- ft[, "Q97.5", ]
wl <- E - 1.959964 * se; wh <- E + 1.959964 * se
cat("cells with lower clamped at 0:", sum(wl < 0), "of", length(wl),
    "; upper clamped at trials:", sum(wh > dm$N), "\n")
cat("cells with lo == 0:", sum(lo == 0), "hi == N:", sum(hi == dm$N), "\n")
k <- which(wl < 0, arr.ind = TRUE)[1:3, , drop = FALSE]
for (r in seq_len(nrow(k))) {
  i <- k[r, 1]; j <- k[r, 2]
  cat(sprintf("row %d cat %d N %d est %.4f se %.4f Wald [%.4f, %.4f] reported [%.4f, %.4f]\n",
              i, j, dm$N[i], E[i, j], se[i, j], wl[i, j], wh[i, j], lo[i, j], hi[i, j]))
}
saveRDS(list(est = ft), file.path("dev/catmulti-rev-out", "frm-mult-vt.rds"))
TRY("fitted linear mult", fitted(fm, scale = "linear"))
TRY("predict mult", predict(fm))
TRY("fitted newdata no trials", fitted(fm, newdata = data.frame(x = 0)))
TRY("fitted newdata trials 5", fitted(fm, newdata = data.frame(x = 0, N = 5)))
TRY("dharma", dharma_residuals(fm))
TRY("ce mult", conditional_effects(fm, "x", categorical = TRUE)[[1]])
ce <- conditional_effects(fm, "x", categorical = TRUE)[[1]]
cat("ce head:\n"); print(head(ce[, intersect(names(ce), c("x", "cats__", "N", "estimate__", "se__"))], 3))
cat("DONE\n")
