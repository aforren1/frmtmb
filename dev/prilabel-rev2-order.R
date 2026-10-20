# Re-check of prilabel fix round 1: brms's group order (frame_re ids)
# against the lane's blocks, re_formula matching, findings 2, 3, 5, 8,
# and frmtmb.sample's unreached-prior note with its absent cases.
# Data dev/dblbar-data.R seed 20261013; extra columns seed 7; gp data
# seed 51 (as dev/prilabel-rev-lscale.R).
.libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
sink(file.path(wt, "dev/prilabel-rev-log/rev2-order.txt"), split = TRUE)
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
dd <- dblbar_data()
set.seed(7)
dd$w1 <- runif(nrow(dd))
dd$w2 <- 1 - dd$w1
dd$ga <- factor(sample(sprintf("a%d", 1:7), nrow(dd), TRUE))
dd$gs <- factor(sample(sprintf("s%d", 1:7), nrow(dd), TRUE))
dd$g1b <- factor(ifelse(as.integer(dd$g1) <= 3L, "lo", "hi"))
dd$yz <- dd$y + rnorm(nrow(dd))
try_ <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
# brms: one row per (id, group, coef, resp, dpar, nlpar) in id order
brms_ids <- function(bf_b) {
  re <- as.data.frame(brms:::frame_re(brms:::brmsterms(bf_b), data = dd))
  paste(re$resp, re$dpar, re$nlpar, re$group, re$coef, sep = "|")
}
# frmtmb: the same, in block order then column order
frm_ids <- function(bf_f) {
  des <- frmtmb:::prior_design(bf_f, dd, NULL, list(),
                               check_response = FALSE)
  out <- character(0)
  for (bk in des$frame$re_blocks) {
    cp <- frmtmb:::block_col_prefix(des$spec, des$frame, bk)
    cf <- frmtmb:::block_col_coefs(bk)
    out <- c(out, paste(cp$resp, cp$dpar, cp$nlpar, bk$group_name, cf,
                        sep = "|"))
  }
  out
}
cases <- list(
  list(y ~ x + (1 | ga) + (1 + x | gr(g, by = hg))),
  list(y ~ x + (1 | ga) + (1 + x || g)),
  list(y ~ x + (1 | ga) + (1 | gr(g, dist = "student"))),
  list(y ~ x + (1 | mm(m1, m2, weights = cbind(w1, w2))) + (1 | g1/g2) +
         (1 | ga), sigma ~ (1 | h) + (1 | g1/g2)),
  list(y ~ x + (1 | p | h) + (1 | ga) + (0 + x | p | h)),
  list(y ~ x + (1 | gs) + (1 | gr(g, by = hg))),
  list(y ~ x + (0 + x | g) + (1 | ga) + (1 | g)),
  list(y ~ x + (1 | ga) + (1 | g:h) + (1 | g))
)
mk <- function(pkg, cs) do.call(getExportedValue(pkg, "bf"), cs)
cat("\n######## group order: brms frame_re ids against frmtmb blocks\n")
n_same <- 0L
for (cs in cases) {
  b <- try_(brms_ids(mk("brms", cs)))
  f <- try_(frm_ids(mk("frmtmb", cs)))
  same <- identical(b, f)
  n_same <- n_same + same
  cat(sprintf("\n%s %s\n", if (same) "SAME" else "DIFF",
              paste(vapply(cs, deparse1, ""), collapse = " ; ")))
  cat(" brms:  ", paste(b, collapse = "  "), "\n")
  if (!same) cat(" frmtmb:", paste(f, collapse = "  "), "\n")
}
cat(sprintf("ORDER SUMMARY %d of %d identical\n", n_same, length(cases)))
mvb <- brms::bf(y ~ x + (1 | ga) + (1 + x | p | g)) +
  brms::bf(yz ~ (1 | p | g) + (1 | gs)) + brms::set_rescor(FALSE)
mvf <- frmtmb::bf(y ~ x + (1 | ga) + (1 + x | p | g)) +
  frmtmb::bf(yz ~ (1 | p | g) + (1 | gs)) + frmtmb::set_rescor(FALSE)
b <- brms_ids(mvb)
f <- try_(frm_ids(mvf))
cat("\nmvbf |p| across responses:", if (identical(b, f)) "SAME" else "DIFF",
    "\n brms:  ", paste(b, collapse = "  "), "\n")
if (!identical(b, f)) cat(" frmtmb:", paste(f, collapse = "  "), "\n")

cat("\n######## re_formula matching on a fit with gr(by), ||, mm, nest\n")
fit <- frm(frmtmb::bf(y ~ x + (1 + x | gr(g, by = hg)) + (1 | ga) +
                        (1 | mm(m1, m2)) + (1 | g1/g2)), data = dd)
nd <- dd[1:4, ]
full <- fitted(fit, newdata = nd)[, 1]
none <- fitted(fit, newdata = nd, re_formula = NA)[, 1]
for (rf in list(~ (1 + x | gr(g, by = hg)), ~ (1 + x | g), ~ (1 | ga),
                ~ (1 | mm(m1, m2)), ~ (1 | g1/g2), ~ (1 | g1:g2),
                ~ (1 | g2:g1), ~ (1 + x | gr(g, by = hg)) + (1 | ga) +
                  (1 | mm(m1, m2)) + (1 | g1/g2))) {
  p <- try_(fitted(fit, newdata = nd, re_formula = rf)[, 1])
  cat(format(deparse1(rf), width = 60), ":",
      if (is.character(p)) substr(p, 1, 110) else
        sprintf("d(full) %.3e d(NA) %.3e", max(abs(p - full)),
                max(abs(p - none))), "\n")
}
fd <- frm(frmtmb::bf(y ~ x + (1 + x || g) + (1 | gr(ga, dist = "student"))),
          data = dd)
for (rf in list(~ (1 + x || g), ~ (1 + x | g), ~ (1 | gr(ga, dist = "student")),
                ~ (1 | ga))) {
  p <- try_(fitted(fd, newdata = nd, re_formula = rf)[, 1])
  cat(format(deparse1(rf), width = 60), ":",
      if (is.character(p)) substr(p, 1, 110) else "runs", "\n")
}

cat("\n######## finding 2: prefix\n")
ent_at <- function(bff, pr) {
  des <- frmtmb:::prior_design(bff, dd, NULL, list(), check_response = FALSE)
  e <- try_(frmtmb:::resolve_priorlist(des, pr)$entries)
  if (is.character(e)) return(substr(e, 1, 120))
  paste(vapply(e, function(z) paste(z$comp, z$idx), ""), collapse = ", ")
}
pr <- function(...) frmtmb::set_prior("normal(0, 0.5)", class = "sd",
                                      group = "g", coef = "x", ...)
f2a <- frmtmb::bf(y ~ x + (1 + x | g), sigma ~ homdiag(1 + x | g))
f2b <- frmtmb::bf(y ~ x + homdiag(1 + x | g), sigma ~ (1 + x | g))
cat("mu us, sigma homdiag, no dpar:   ", ent_at(f2a, pr()), "\n")
cat("mu us, sigma homdiag, dpar sigma:", ent_at(f2a, pr(dpar = "sigma")), "\n")
cat("mu homdiag, sigma us, no dpar:   ", ent_at(f2b, pr()), "\n")
cat("mu homdiag, sigma us, dpar sigma:", ent_at(f2b, pr(dpar = "sigma")), "\n")
cat("mv: y homdiag, yz us, resp yz:   ",
    ent_at(frmtmb::bf(y ~ homdiag(1 + x | g)) + frmtmb::bf(yz ~ (1 + x | g)) +
             frmtmb::set_rescor(FALSE), pr(resp = "yz")), "\n")
cat("mv: y homdiag, yz us, resp y:    ",
    ent_at(frmtmb::bf(y ~ homdiag(1 + x | g)) + frmtmb::bf(yz ~ (1 + x | g)) +
             frmtmb::set_rescor(FALSE), pr(resp = "y")), "\n")

cat("\n######## finding 3: sds_ after cor_\n")
set.seed(61)
d3 <- data.frame(x = runif(240), g = factor(sample(1:10, 240, TRUE)),
                 h = factor(sample(1:6, 240, TRUE)))
d3$y <- sin(4 * d3$x) + rnorm(10, 0, 0.6)[d3$g] +
  rnorm(10, 0, 0.4)[d3$g] * d3$x + rnorm(6, 0, 0.4)[d3$h] + rnorm(240, 0, 0.3)
for (f in list(y ~ s(x) + (1 + x | g), y ~ (1 | h) + s(x) + (1 + x | g) +
                 (0 + x | h))) {
  v <- variables(frm(frmtmb::bf(f), data = d3))
  cat(deparse1(f), ":\n ", v[!grepl("^b_|^r_|^s_|^zs_|^bs_", v)], "\n")
}
cat("brms (round 0, order-brms.txt): sd_g__Intercept sd_g__x sd_h__Intercept",
    "sd_h__x cor_g__Intercept__x sds_sx_1 sigma\n")

cat("\n######## finding 5: levels with ':'\n")
set.seed(9)
dd$k1 <- factor(sample(c("c:d", "e"), nrow(dd), TRUE))
dd$k2 <- factor(sample(c("f", "h"), nrow(dd), TRUE))
dd$k3 <- factor(sample(c("u:v", "w"), nrow(dd), TRUE))
f5 <- frm(frmtmb::bf(y ~ x + (1 | k1/k2) + (1 | k2:k1) + (1 | k1:k2:k3)),
          data = dd)
print(lapply(ranef(f5), function(a) dimnames(a)[[1L]]))
print(lapply(coef(f5)[[1L]], function(a) rownames(a)))
nd5 <- dd[1:2, ]
nd5$k2 <- factor("zz")
cat(try_(fitted(f5, newdata = nd5)), "\n")
b5 <- brms:::combine_groups(dd, "k1:k2", "k2:k1", "k1:k2:k3")
cat("brms combine_groups:", sort(unique(b5[["k1:k2"]])), "|",
    sort(unique(b5[["k2:k1"]])), "|", sort(unique(b5[["k1:k2:k3"]])), "\n")
pf <- fitted(f5, newdata = dd[1:5, ])[, 1]
cat(sprintf("newdata at old levels against fitted: %.3e\n",
            max(abs(pf - fitted(f5)[1:5, 1]))))

cat("\n######## finding 8: validate_prior sd order, mv\n")
pb <- brms::set_prior("normal(0, 0.5)", class = "sd", group = "ga",
                      resp = "y")
vb <- as.data.frame(brms::validate_prior(pb, mvb, data = dd))
vf <- as.data.frame(frmtmb::validate_prior(
  frmtmb::set_prior("normal(0, 0.5)", class = "sd", group = "ga", resp = "y"),
  mvf, data = dd))
k <- c("class", "coef", "group", "resp")
sb <- vb[vb$class == "sd", k]
sf <- vf[vf$class == "sd", k]
rownames(sb) <- rownames(sf) <- NULL
cat("sd rows identical:", identical(sb, sf), "\n")
if (!identical(sb, sf)) { print(sb); print(sf) }

cat("\n######## finding 7: frm_sample note for a prior that reaches nothing\n")
set.seed(51)
n <- 80
dg <- data.frame(x = runif(n), z = runif(n), g = factor(sample(1:8, n, TRUE)))
dg$y <- sin(4 * dg$x) + rnorm(8, 0, 0.5)[dg$g] + rnorm(n, 0, 0.3)
note_for <- function(form, pr, data = dg) {
  uf <- frm(frmtmb::bf(form), data = data, dry_run = "objective")
  msgs <- character(0)
  r <- withCallingHandlers(
    try_(frmtmb.sample:::sample_resolve_priors(uf, pr)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  hit <- grep("applies to no parameter", unlist(strsplit(msgs, "\n")),
              value = TRUE)
  if (is.character(r)) return(substr(r, 1, 120))
  if (length(hit)) paste("NOTE:", trimws(hit)) else "no note"
}
sp <- frmtmb::set_prior
gpf <- y ~ gp(x) + (1 | g)
cat("class lscale (presence):     ", note_for(gpf, sp("inv_gamma(2, 1)",
                                                      class = "lscale")),
    "\n")
cat("lscale coef gpx:             ", note_for(gpf, sp("inv_gamma(2, 1)",
                                                      class = "lscale",
                                                      coef = "gpx")), "\n")
cat("sd group g:                  ", note_for(gpf, sp("normal(0, 1)",
                                                      class = "sd",
                                                      group = "g")), "\n")
cat("sd class-wide:               ", note_for(gpf, sp("normal(0, 1)",
                                                      class = "sd")), "\n")
cat("sdgp class-wide:             ", note_for(gpf, sp("normal(0, 1)",
                                                      class = "sdgp")), "\n")
cat("sigma natural:               ", note_for(gpf, sp("student_t(3, 0, 1)",
                                                      class = "sigma")), "\n")
cat("b on x (flat default):       ", note_for(y ~ x + (1 | g),
                                              sp("normal(0, 1)",
                                                 class = "b")), "\n")
cat("cor lkj(2):                  ", note_for(y ~ x + (1 + x | g),
                                              sp("lkj(2)", class = "cor")),
    "\n")
cat("sd with lb:                  ", note_for(y ~ x + (1 | g),
                                              sp("normal(0, 1)", class = "sd",
                                                 lb = 0.01)), "\n")
cat("sd class + coef on all coefs:", note_for(y ~ x + (1 | g),
                                              sp("normal(0, 1)", class = "sd") +
                                                sp("normal(0, 2)", class = "sd",
                                                   group = "g",
                                                   coef = "Intercept")), "\n")
cat("two gp terms, class lscale:  ", note_for(y ~ gp(x) + gp(z),
                                              sp("inv_gamma(2, 1)",
                                                 class = "lscale")), "\n")
sink()
