# Reviewer: try to falsify the cmc claims against brms::standata() and
# default_prior() on the lane's build.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) structure(conditionMessage(e), class = "ERR"))
}
set.seed(21)
n <- 240
d <- data.frame(f = factor(sample(c("a", "b", "c"), n, TRUE)),
                f2 = factor(sample(c("p", "q"), n, TRUE)),
                x = rnorm(n), h = factor(rep(1:12, each = 20)))
d$y <- rnorm(n, as.numeric(d$f) + d$x, 0.7) +
  rnorm(12, 0, 0.5)[as.integer(d$h)]
d$y2 <- d$y + rnorm(n)
d$ord <- cut(d$y, quantile(d$y, c(0, .3, .6, 1)), include.lowest = TRUE,
             labels = FALSE)

strip <- function(X) { X <- as.matrix(X); dimnames(X) <- NULL
  attr(X, "assign") <- NULL; attr(X, "contrasts") <- NULL; X }
key <- function(fr, dp = "mu", resp = "y") paste0(resp, ".", dp)

cmp_X <- function(label, ffrm, fbrms, dp = "mu", sdname = "X",
                  resp = "y", family = NULL, bfam = brms::brmsfamily("gaussian")) {
  fr <- q(frm(ffrm, data = d, family = family, dry_run = "frame"))
  sd <- q(brms::standata(fbrms, d, family = bfam))
  if (inherits(fr, "ERR") || inherits(sd, "ERR")) {
    cat(sprintf("%-38s frm: %s | brms: %s\n", label,
                if (inherits(fr, "ERR")) paste("ERR", substr(fr, 1, 60))
                else "ok",
                if (inherits(sd, "ERR")) paste("ERR", substr(sd, 1, 60))
                else "ok"))
    return(invisible())
  }
  lp <- fr$linpreds[[paste0(resp, ".", dp)]]
  if (is.null(lp)) lp <- fr$linpreds[[1]]
  Xf <- lp$X[, seq_len(lp$n_param_cols %||% ncol(lp$X)), drop = FALSE]
  Xb <- sd[[sdname]]
  same_cols <- identical(colnames(Xf), colnames(Xb))
  same_val <- same_cols && isTRUE(all.equal(strip(Xf), strip(Xb)))
  cat(sprintf("%-38s cols frm [%s] brms [%s] identical-values %s\n", label,
              paste(colnames(Xf), collapse = ","),
              paste(colnames(Xb), collapse = ","), same_val))
  invisible(list(fr = fr, sd = sd))
}
B <- brms::bf
cat("== population-level\n")
for (cm in c(TRUE, FALSE)) {
  cat("-- cmc =", cm, "\n")
  cmp_X("0 + f", bf(y ~ 0 + f, cmc = cm), B(y ~ 0 + f, cmc = cm))
  cmp_X("0 + f:x", bf(y ~ 0 + f:x, cmc = cm), B(y ~ 0 + f:x, cmc = cm))
  cmp_X("0 + f + f2", bf(y ~ 0 + f + f2, cmc = cm),
        B(y ~ 0 + f + f2, cmc = cm))
  cmp_X("0 + f * x", bf(y ~ 0 + f * x, cmc = cm), B(y ~ 0 + f * x, cmc = cm))
  cmp_X("0 + x + f", bf(y ~ 0 + x + f, cmc = cm), B(y ~ 0 + x + f, cmc = cm))
  cmp_X("0 + Intercept + f", bf(y ~ 0 + Intercept + f, cmc = cm),
        B(y ~ 0 + Intercept + f, cmc = cm))
  cmp_X("-1 + f", bf(y ~ -1 + f, cmc = cm), B(y ~ -1 + f, cmc = cm))
  cmp_X("sigma: lf(0 + f)", bf(y ~ 1) + lf(sigma ~ 0 + f, cmc = cm),
        B(y ~ 1) + brms::lf(sigma ~ 0 + f, cmc = cm), dp = "sigma",
        sdname = "X_sigma")
  cmp_X("bf(cmc) does not reach sigma",
        bf(y ~ 0 + f, sigma ~ 0 + f, cmc = cm),
        B(y ~ 0 + f, sigma ~ 0 + f, cmc = cm), dp = "sigma",
        sdname = "X_sigma")
  cmp_X("cumulative ord ~ 0 + f", bf(ord ~ 0 + f, cmc = cm),
        B(ord ~ 0 + f, cmc = cm), resp = "ord", family = cumulative(),
        bfam = brms::cumulative())
  cmp_X("nlpar lf(a ~ 0 + f)",
        bf(y ~ a + b * x, b ~ 1, nl = TRUE) + lf(a ~ 0 + f, cmc = cm),
        B(y ~ a + b * x, b ~ 1, nl = TRUE) + brms::lf(a ~ 0 + f, cmc = cm),
        dp = "a", sdname = "X_a")
}

cat("== group-level: Z of brms vs frame's component cnms\n")
cmpZ <- function(label, ffrm, fbrms) {
  fr <- q(frm(ffrm, data = d, dry_run = "frame"))
  sd <- q(brms::standata(fbrms, d))
  if (inherits(fr, "ERR") || inherits(sd, "ERR")) {
    cat(label, "frm", if (inherits(fr, "ERR")) fr else "ok", "brms",
        if (inherits(sd, "ERR")) sd else "ok", "\n")
    return(invisible())
  }
  zb <- names(sd)[grepl("^Z_", names(sd))]
  cn <- unlist(lapply(fr$re_blocks, function(b)
    lapply(b$components, `[[`, "cnms")))
  # brms's Z_k_j for one term, column j
  lp <- fr$linpreds[["y.mu"]]
  Z <- as.matrix(lp$Z)
  nc <- length(cn)
  ok <- nc == length(zb) &&
    all(vapply(seq_len(nc), function(j) isTRUE(all.equal(
      unname(rowSums(Z[, seq(j, ncol(Z), by = nc), drop = FALSE])),
      as.numeric(sd[[zb[j]]]))), NA))
  cat(sprintf("%-34s frm cnms [%s]  brms Z [%s]  values %s\n", label,
              paste(cn, collapse = ","), paste(zb, collapse = ","), ok))
}
for (cm in c(TRUE, FALSE)) {
  cat("-- cmc =", cm, "\n")
  cmpZ("(0 + f | h)", bf(y ~ x + (0 + f | h), cmc = cm),
       B(y ~ x + (0 + f | h), cmc = cm))
  cmpZ("0 + x + (0 + f | h)", bf(y ~ 0 + x + (0 + f | h), cmc = cm),
       B(y ~ 0 + x + (0 + f | h), cmc = cm))
  cmpZ("(0 + f + x | h)", bf(y ~ x + (0 + f + x | h), cmc = cm),
       B(y ~ x + (0 + f + x | h), cmc = cm))
  cmpZ("(0 + x | h)", bf(y ~ x + (0 + x | h), cmc = cm),
       B(y ~ x + (0 + x | h), cmc = cm))
  cmpZ("(0 + f:x | h)", bf(y ~ x + (0 + f:x | h), cmc = cm),
       B(y ~ x + (0 + f:x | h), cmc = cm))
  cmpZ("(0 + f || h)", bf(y ~ x + (0 + f || h), cmc = cm),
       B(y ~ x + (0 + f || h), cmc = cm))
  cmpZ("(f | h)", bf(y ~ x + (f | h), cmc = cm),
       B(y ~ x + (f | h), cmc = cm))
}

cat("== default_prior rows, cmc = FALSE, against brms\n")
dp_f <- q(default_prior(bf(y ~ 0 + f + (0 + f | h), cmc = FALSE), data = d))
dp_b <- q(brms::default_prior(B(y ~ 0 + f + (0 + f | h), cmc = FALSE),
                               data = d))
print(dp_f[, c("class", "coef", "group")])
print(as.data.frame(dp_b)[, c("prior", "class", "coef", "group")])

cat("== a fit: names, prediction at newdata, unseen levels\n")
fit <- q(frm(bf(y ~ 0 + f + (0 + f | h), cmc = FALSE), data = d))
print(fixef(fit)); print(variables(fit)); str(ranef(fit))

nd <- d[1:6, ]
lin <- q(frm_linpred(fit, newdata = nd))
cat("linpred at newdata == in-sample:",
    isTRUE(all.equal(unname(lin), unname(frm_linpred(fit)[1:6]))), "\n")
# hand-computed eta: X beta + Z b
be <- fixef(fit)[, "Estimate"]
re <- ranef(fit)$h
eta <- be["fb"] * (nd$f == "b") + be["fc"] * (nd$f == "c") +
  re[as.character(nd$h), "fb"] * (nd$f == "b") +
  re[as.character(nd$h), "fc"] * (nd$f == "c")
cat("hand eta vs frm_linpred max abs diff / sd(y):",
    max(abs(unname(eta) - unname(lin))) / sd(d$y), "\n")
nd_new <- nd; nd_new$h <- factor(c("new1", "new1", "new2", "new2", "new3",
                                   "new3"))
p1 <- q(predict(fit, newdata = nd_new, allow_new_levels = TRUE))
cat("new group level:", if (inherits(p1, "ERR")) p1 else "ok", "\n")
p2 <- q(fitted(fit, newdata = nd_new, re_formula = NA))
cat("re_formula = NA:", if (inherits(p2, "ERR")) p2 else "ok", "\n")
nd_f <- nd; nd_f$f <- factor(c("a", "b", "c", "d", "a", "b"))
p3 <- q(fitted(fit, newdata = nd_f))
cat("unseen fixed factor level d:", if (inherits(p3, "ERR")) p3 else
  "ok (no refusal)", "\n")
p4 <- q(fitted(fit, newdata = nd, re_formula = ~ (0 + f | h)))
p0 <- q(fitted(fit, newdata = nd))
cat("re_formula = ~ (0 + f | h):",
    if (inherits(p4, "ERR")) p4 else
      isTRUE(all.equal(p4, p0)), "\n")
nd_one <- nd[nd$f == "a", ][1:2, ]
nd_one$f <- factor(nd_one$f)   # a newdata factor with one level
p5 <- q(fitted(fit, newdata = nd_one))
cat("newdata factor with only level a:",
    if (inherits(p5, "ERR")) p5 else "ok", "\n")
s1 <- q(simulate(fit, nsim = 2, seed = 1, newdata = nd))
cat("simulate(newdata):", if (inherits(s1, "ERR")) s1 else "ok", "\n")
ce <- q(conditional_effects(fit, "f"))
cat("conditional_effects:", if (inherits(ce, "ERR")) ce else "ok", "\n")
em <- q(emmeans::emmeans(fit, ~ f))
cat("emmeans:", if (inherits(em, "ERR")) em else "ok", "\n")
if (!inherits(em, "ERR")) print(em)
mm <- q(model.matrix(fit))
cat("model.matrix cols:", if (inherits(mm, "ERR")) mm else colnames(mm),
    "\n")

cat("== update\n")
f0 <- q(frm(y ~ 0 + f, data = d))
u <- q(update(f0, bf(~ ., cmc = FALSE)))
cat("update(f0, bf(~ ., cmc = FALSE)) fixef:",
    if (inherits(u, "ERR")) u else rownames(fixef(u)), "\n")
u2 <- q(update(fit, newdata = d[-(1:5), ]))
cat("update newdata fixef:", if (inherits(u2, "ERR")) u2 else
  rownames(fixef(u2)), "\n")
u3 <- q(update(fit, . ~ . + x))
cat("update . ~ . + x fixef:", if (inherits(u3, "ERR")) u3 else
  rownames(fixef(u3)), "\n")

cat("== multi-membership and gr(by =)\n")
d$h2 <- factor(sample(levels(d$h), n, TRUE))
cmp_mm <- q(frm(bf(y ~ x + (0 + f | mm(h, h2)), cmc = FALSE), data = d,
                dry_run = "frame"))
if (!inherits(cmp_mm, "ERR")) {
  cat("mm cnms:", unlist(lapply(cmp_mm$re_blocks, function(b)
    lapply(b$components, `[[`, "cnms"))), "\n")
} else cat("mm:", cmp_mm, "\n")
sdmm <- q(brms::standata(B(y ~ x + (0 + f | mm(h, h2)), cmc = FALSE), d))
cat("brms mm Z:", if (inherits(sdmm, "ERR")) sdmm else
  grep("^Z_", names(sdmm), value = TRUE), "\n")
d$hb <- factor(ifelse(as.integer(d$h) %% 2 == 0, "e", "o"))
cby <- q(frm(bf(y ~ x + (0 + f | gr(h, by = hb)), cmc = FALSE), data = d,
             dry_run = "frame"))
cat("gr(by) cnms:", if (inherits(cby, "ERR")) cby else
  unlist(lapply(cby$re_blocks, function(b)
    lapply(b$components, `[[`, "cnms"))), "\n")
sdby <- q(brms::standata(B(y ~ x + (0 + f | gr(h, by = hb)), cmc = FALSE),
                         d))
cat("brms gr(by) Z:", if (inherits(sdby, "ERR")) sdby else
  grep("^Z_", names(sdby), value = TRUE), "\n")
cat("DONE\n")
