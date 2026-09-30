# Reviewer re-check (lane ceplot, punch round 1): "old_levels" on shapes
# the first round did not try, against brms 2.23.0 at seeds 1 to 12:
#   (e) (1 | g) + (1 | g:h), g unseen (so g:h unseen too)
#   (f) mm(g1, g2): one member unseen; both at one unseen value; two
#       different unseen values
#   (g) a multivariate model whose two responses both read (1 | g)
# and, for frmtmb, that predict() makes the same choice as fitted() at
# the same seed. The choice is DECODED from the answer at x = 0: the
# prediction minus the fixed part is a sum of seen levels' effects, and
# the (unique) combination that reproduces it to 1e-9 relative is the
# choice. brms: fixed_param, 40 chains of one iteration at random inits.
# Data seed 3.
#   Rscript dev/ceplot-rev-oldlevels3.R > dev/ceplot-rev-log/p1/oldlevels3.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
})
options(mc.cores = 1)
say <- function(...) cat(sprintf(...), "\n", sep = "")
qb <- function(e) suppressMessages(suppressWarnings(e))
bfix <- function(form, data) {
  qb(brm(form, data = data, family = stats::gaussian(),
         algorithm = "fixed_param", chains = 40, iter = 1, warmup = 0,
         refresh = 0, seed = 1, silent = 2))
}
# decode: which candidate combination (rows of `cand`, each a list of
# per-term effect vectors) reproduces `resid`
decode <- function(resid, cand) {
  err <- vapply(cand$eff, function(e) max(abs(e - resid)), 1)
  k <- which(err <= 1e-9 * max(1, max(abs(resid))))
  if (length(k) == 1L) cand$lab[k] else if (length(k)) "AMBIG" else "none"
}
set.seed(3)
ng <- 12L
d <- data.frame(x = rnorm(360), g = factor(rep(seq_len(ng), 30)),
                h = factor(rep(1:3, each = 120)))
d$g1 <- factor(sample(1:ng, 360, TRUE), levels = 1:ng)
d$g2 <- factor(sample(1:ng, 360, TRUE), levels = 1:ng)
u <- rnorm(ng, 0, 2)
d$y <- rnorm(360, 1 + 0.5 * d$x + u[d$g] +
               rnorm(36)[as.integer(interaction(d$g, d$h))])
d$ym <- rnorm(360, 1 + 0.5 * d$x + 0.5 * (u[d$g1] + u[d$g2]))
d$y2 <- rnorm(360, -1 + 0.3 * d$x + u[d$g])

## (e) nested
fe <- frm(bf(y ~ x + (1 | g) + (1 | g:h)), family = gaussian(), data = d)
be <- bfix(y ~ x + (1 | g) + (1 | g:h), d)
re <- ranef(fe)
rg <- setNames(re$g[, "Intercept"], rownames(re$g))
rgh <- setNames(re$`g:h`[, "Intercept"], rownames(re$`g:h`))
cf <- expand.grid(a = names(rg), b = names(rgh), stringsAsFactors = FALSE)
cand_f <- list(eff = lapply(seq_len(nrow(cf)), function(i) {
  rg[[cf$a[i]]] + rgh[[cf$b[i]]]
}), lab = paste0("g=", cf$a, ",g:h=", cf$b))
bm <- as_draws_matrix(be)
bg <- bm[, grep("^r_g\\[", colnames(bm)), drop = FALSE]
bgh <- bm[, grep("^r_g:h\\[", colnames(bm)), drop = FALSE]
lg <- sub("^r_g\\[(.*),Intercept\\]$", "\\1", colnames(bg))
lgh <- sub("^r_g:h\\[(.*),Intercept\\]$", "\\1", colnames(bgh))
cb <- expand.grid(a = seq_along(lg), b = seq_along(lgh))
cand_b <- list(eff = lapply(seq_len(nrow(cb)), function(i) {
  as.numeric(bg[, cb$a[i]] + bgh[, cb$b[i]])
}), lab = paste0("g=", lg[cb$a], ",g:h=", sub("_", ":", lgh[cb$b])))
nd <- data.frame(x = 0, g = factor("n1"), h = factor("1", levels = 1:3))
fix_f <- fixef(fe)["Intercept", "Estimate"]
res_e <- t(vapply(1:12, function(s) {
  set.seed(s)
  v <- fitted(fe, newdata = nd, allow_new_levels = TRUE,
              sample_new_levels = "old_levels")[1, 1]
  set.seed(s)
  vb <- fitted(be, newdata = nd, allow_new_levels = TRUE,
               sample_new_levels = "old_levels", summary = FALSE)[, 1]
  c(decode(v - fix_f, cand_f),
    decode(vb - as.numeric(bm[, "b_Intercept"]), cand_b))
}, character(2)))
say("(e) frmtmb: %s", paste(res_e[, 1], collapse = " | "))
say("(e) brms:   %s", paste(res_e[, 2], collapse = " | "))
say("(e) agree at %d of 12", sum(res_e[, 1] == res_e[, 2]))
# predict() makes fitted()'s choice at the same seed
same_pred <- vapply(1:12, function(s) {
  set.seed(s)
  a <- frmtmb:::predict_new_level_spec(fe, fe$spec$responses[[1L]], nd,
                                       NULL, TRUE, "old_levels")
  set.seed(s)
  o <- frmtmb:::fitted_old_levels(fe, names(fe$spec$responses)[1L], nd,
                                  NULL)
  identical(unlist(attr(a, "old_pick")), unlist(o[["new_level_pick"]]))
}, NA)
say("(e) predict()'s choice equals fitted()'s at %d of 12 seeds",
    sum(same_pred))

## (f) mm
fm <- frm(bf(ym ~ x + (1 | mm(g1, g2))), family = gaussian(), data = d)
bmm <- bfix(ym ~ x + (1 | mm(g1, g2)), d)
rf <- ranef(fm)[[1]][, "Intercept"]
names(rf) <- rownames(ranef(fm)[[1]])
bmd <- as_draws_matrix(bmm)
br <- bmd[, grep("^r_mm", colnames(bmd)), drop = FALSE]
lb <- sub("^r_[^[]*\\[(.*),Intercept\\]$", "\\1", colnames(br))
mm_case <- function(g1, g2, seen2) {
  nd <- data.frame(x = 0, g1 = factor(g1), g2 = factor(g2))
  # candidates: members' levels k1 (and k2 when member 2 is unseen)
  ks <- names(rf)
  pairs <- if (seen2) expand.grid(a = ks, b = g2, stringsAsFactors = FALSE)
           else expand.grid(a = ks, b = ks, stringsAsFactors = FALSE)
  cf <- list(eff = lapply(seq_len(nrow(pairs)), function(i) {
    0.5 * rf[[pairs$a[i]]] + 0.5 * rf[[pairs$b[i]]]
  }), lab = paste0(pairs$a, "+", pairs$b))
  cbb <- list(eff = lapply(seq_len(nrow(pairs)), function(i) {
    as.numeric(0.5 * br[, match(pairs$a[i], lb)] +
                 0.5 * br[, match(pairs$b[i], lb)])
  }), lab = cf$lab)
  t(vapply(1:12, function(s) {
    set.seed(s)
    v <- tryCatch(fitted(fm, newdata = nd, allow_new_levels = TRUE,
                         sample_new_levels = "old_levels")[1, 1],
                  error = function(e) NA_real_)
    set.seed(s)
    vb <- tryCatch(fitted(bmm, newdata = nd, allow_new_levels = TRUE,
                          sample_new_levels = "old_levels",
                          summary = FALSE)[, 1], error = function(e) NULL)
    c(if (is.na(v)) "ERR" else
        decode(v - fixef(fm)["Intercept", "Estimate"], cf),
      if (is.null(vb)) "ERR" else
        decode(vb - as.numeric(bmd[, "b_Intercept"]), cbb))
  }, character(2)))
}
for (cs in list(list("n1", "2", TRUE, "one member unseen"),
                list("n1", "n1", FALSE, "both at one unseen value"),
                list("n1", "n2", FALSE, "two unseen values"))) {
  r <- mm_case(cs[[1]], cs[[2]], cs[[3]])
  say("(f) %s: frmtmb %s", cs[[4]], paste(r[, 1], collapse = " "))
  say("(f) %s: brms   %s", cs[[4]], paste(r[, 2], collapse = " "))
  say("(f) %s: agree at %d of 12", cs[[4]], sum(r[, 1] == r[, 2]))
}

## (g) multivariate, (1 | g) in both responses
fv <- frm(bf(y ~ x + (1 | g)) + bf(y2 ~ x + (1 | g)), family = gaussian(),
          data = d)
bv <- qb(brm(brms::bf(y ~ x + (1 | g)) + brms::bf(y2 ~ x + (1 | g)) +
               brms::set_rescor(FALSE), data = d,
             family = stats::gaussian(), algorithm = "fixed_param",
             chains = 40, iter = 1, warmup = 0, refresh = 0, seed = 1,
             silent = 2))
rv <- ranef(fv)
nd <- data.frame(x = 0, g = factor("n1"))
pick_f <- function(v, r, int) {
  lv <- rownames(r)
  e <- v - int
  k <- which(abs(r[, "Intercept"] - e) <= 1e-9 * max(1, abs(e)))
  if (length(k) == 1L) lv[k] else "none"
}
bvd <- as_draws_matrix(bv)
pick_b <- function(v, resp) {
  cols <- grep(sprintf("^r_g__%s\\[", resp), colnames(bvd), value = TRUE)
  e <- v - as.numeric(bvd[, sprintf("b_%s_Intercept", resp)])
  k <- which(vapply(cols, function(cn) max(abs(bvd[, cn] - e)) < 1e-9, NA))
  if (length(k) == 1L) sub("^.*\\[(.*),Intercept\\]$", "\\1", cols[k]) else
    "none"
}
say("(g) frmtmb ranef names: %s", paste(names(rv), collapse = ","))
res_g <- t(vapply(1:12, function(s) {
  set.seed(s)
  f <- fitted(fv, newdata = nd, allow_new_levels = TRUE,
              sample_new_levels = "old_levels")
  set.seed(s)
  b <- fitted(bv, newdata = nd, allow_new_levels = TRUE,
              sample_new_levels = "old_levels", summary = FALSE)
  fe1 <- f[1, "Estimate", 1]
  fe2 <- f[1, "Estimate", 2]
  c(pick_f(fe1, rv[[1]], fixef(fv)["y_Intercept", "Estimate"]),
    pick_f(fe2, rv[[2]], fixef(fv)["y2_Intercept", "Estimate"]),
    pick_b(b[, 1, 1], "y"), pick_b(b[, 1, 2], "y2"))
}, character(4)))
say("(g) frmtmb y/y2: %s", paste(res_g[, 1], res_g[, 2], sep = "/",
                                 collapse = " "))
say("(g) brms   y/y2: %s", paste(res_g[, 3], res_g[, 4], sep = "/",
                                 collapse = " "))
say("(g) frmtmb one group for both responses at %d of 12; brms at %d",
    sum(res_g[, 1] == res_g[, 2]), sum(res_g[, 3] == res_g[, 4]))
say("done")
