# Reviewer check (lane ceplot): fitted(sample_new_levels = "old_levels")
# against brms 2.23.0 on shapes the worker did not try: which seen level
# each unseen level reads, at seeds 1 to 12, in
#   (a) two different unseen levels in one newdata;
#   (b) (1 | g) in mu AND in sigma (one grouping factor, two blocks),
#       read through dpar = "mu" and dpar = "sigma";
#   (c) gr(g, by = f): the choice must stay within the row's by-level;
#   (d) (1 | g) + (1 | h), both unseen.
# brms: algorithm = "fixed_param", 40 chains of one iteration at random
# inits, so its seen levels differ; its choice is read by matching its
# own fitted() at each seen level (all draws identical). frmtmb: the
# estimate identical() to fitted() at a seen level.
# Data seed 3 (a, b, d) and 4 (c).
#   Rscript dev/ceplot-rev-oldlevels.R > dev/ceplot-rev-log/oldlevels.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
})
options(mc.cores = 1)
cat("frmtmb from", find.package("frmtmb"), "\n")
say <- function(...) cat(sprintf(...), "\n", sep = "")
qb <- function(e) suppressMessages(suppressWarnings(e))
bfix <- function(form, data, family = stats::gaussian()) {
  qb(brm(form, data = data, family = family, algorithm = "fixed_param",
         chains = 40, iter = 1, warmup = 0, refresh = 0, seed = 1,
         silent = 2))
}
# which seen level (row of `levs`) a call's rows read, by exact match
which_level <- function(val, cand) {
  hit <- which(vapply(cand, function(a) isTRUE(all.equal(a, val,
                                                         tolerance = 0)),
                      NA))
  if (length(hit) == 1L) hit else NA_integer_
}
set.seed(3)
ng <- 12L
dg <- data.frame(x = rnorm(240), g = factor(rep(seq_len(ng), 20)),
                 h = factor(rep(1:8, 30)))
dg$y <- rnorm(240, 1 + 0.5 * dg$x + rnorm(ng, 0, 2)[dg$g] +
                rnorm(8)[dg$h], exp(0.2 + rnorm(ng, 0, 0.5)[dg$g]))
dg$yl <- exp(dg$y / 4)

## (a) two unseen levels in one newdata
fa <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dg)
ba <- bfix(y ~ x + (1 | g), dg)
nd <- data.frame(x = 0, g = factor(c("n1", "n2", "n1")))
seen_f <- lapply(seq_len(ng), function(k) {
  unname(fitted(fa, newdata = data.frame(x = 0, g = factor(k, levels = 1:ng)))[1, 1:2])
})
seen_b <- lapply(seq_len(ng), function(k) {
  fitted(ba, newdata = data.frame(x = 0, g = factor(k, levels = 1:ng)),
         summary = FALSE)[, 1]
})
pf <- pb <- matrix(NA_integer_, 12, 3)
for (s in 1:12) {
  set.seed(s)
  f <- fitted(fa, newdata = nd, allow_new_levels = TRUE,
              sample_new_levels = "old_levels")
  pf[s, ] <- vapply(1:3, function(i) which_level(unname(f[i, 1:2]), seen_f),
                    1L)
  set.seed(s)
  b <- fitted(ba, newdata = nd, allow_new_levels = TRUE,
              sample_new_levels = "old_levels", summary = FALSE)
  pb[s, ] <- vapply(1:3, function(i) which_level(b[, i], seen_b), 1L)
}
say("(a) frmtmb rows n1,n2,n1 per seed: %s",
    paste(apply(pf, 1, paste, collapse = "/"), collapse = " "))
say("(a) brms   rows n1,n2,n1 per seed: %s",
    paste(apply(pb, 1, paste, collapse = "/"), collapse = " "))
say("(a) seeds where every row agrees: %d of 12", sum(rowSums(pf == pb) == 3,
                                                     na.rm = TRUE))

## (b) (1 | g) in mu and in sigma, lognormal so the mean reads both
fb <- frm(bf(yl ~ x + (1 | g), sigma ~ (1 | g)), family = lognormal(),
          data = dg)
bb <- bfix(brms::bf(yl ~ x + (1 | g), sigma ~ (1 | g)), dg, brms::lognormal())
nd1 <- data.frame(x = 0, g = factor("n1"))
one <- function(k) data.frame(x = 0, g = factor(k, levels = 1:ng))
pick <- function(obj, s, dp, brms_side) {
  set.seed(s)
  if (brms_side) {
    v <- fitted(obj, newdata = nd1, allow_new_levels = TRUE,
                sample_new_levels = "old_levels", dpar = dp,
                summary = FALSE)[, 1]
    which_level(v, lapply(seq_len(ng), function(k) {
      fitted(obj, newdata = one(k), dpar = dp, summary = FALSE)[, 1]
    }))
  } else {
    v <- unname(fitted(obj, newdata = nd1, allow_new_levels = TRUE,
                       sample_new_levels = "old_levels", dpar = dp)[, 1])
    which_level(v, lapply(seq_len(ng), function(k) {
      unname(fitted(obj, newdata = one(k), dpar = dp)[, 1])
    }))
  }
}
res <- t(vapply(1:12, function(s) {
  c(pick(fb, s, "mu", FALSE), pick(fb, s, "sigma", FALSE),
    pick(bb, s, "mu", TRUE), pick(bb, s, "sigma", TRUE))
}, integer(4)))
say("(b) frmtmb mu/sigma level per seed: %s",
    paste(apply(res[, 1:2], 1, paste, collapse = "/"), collapse = " "))
say("(b) brms   mu/sigma level per seed: %s",
    paste(apply(res[, 3:4], 1, paste, collapse = "/"), collapse = " "))
say("(b) frmtmb reads ONE seen group for both at %d of 12 seeds; brms at %d",
    sum(res[, 1] == res[, 2], na.rm = TRUE),
    sum(res[, 3] == res[, 4], na.rm = TRUE))
# the response-scale mean: exp(mu + sigma^2 / 2) at a mixed pair is no
# seen group's mean
set.seed(1)
fm <- fitted(fb, newdata = nd1, allow_new_levels = TRUE,
             sample_new_levels = "old_levels")[, 1]
means <- vapply(seq_len(ng), function(k) fitted(fb, newdata = one(k))[, 1], 1)
say("(b) seed 1 response mean %.5f; equals a seen group's mean: %s",
    fm, any(means == fm))

## (c) gr(g, by = f): f constant within g
set.seed(4)
dc <- data.frame(x = rnorm(240), g = factor(rep(seq_len(ng), 20)))
dc$f <- factor(ifelse(as.integer(dc$g) <= 6, "a", "b"))
dc$y <- rnorm(240, 1 + 0.5 * dc$x + rnorm(ng, 0, 2)[dc$g])
fc <- frm(bf(y ~ x + (1 | gr(g, by = f))), family = gaussian(), data = dc)
bc <- bfix(y ~ x + (1 | gr(g, by = f)), dc)
ndc <- data.frame(x = 0, g = factor("n1"), f = factor("b", levels = c("a", "b")))
onec <- function(k) data.frame(x = 0, g = factor(k, levels = 1:ng),
                               f = factor(ifelse(k <= 6, "a", "b"),
                                          levels = c("a", "b")))
cf <- vapply(1:12, function(s) {
  set.seed(s)
  v <- unname(fitted(fc, newdata = ndc, allow_new_levels = TRUE,
                     sample_new_levels = "old_levels")[, 1:2, drop = FALSE])
  which_level(v, lapply(seq_len(ng), function(k) {
    unname(fitted(fc, newdata = onec(k))[, 1:2, drop = FALSE])
  }))
}, 1L)
cb <- vapply(1:12, function(s) {
  set.seed(s)
  v <- fitted(bc, newdata = ndc, allow_new_levels = TRUE,
              sample_new_levels = "old_levels", summary = FALSE)[, 1]
  which_level(v, lapply(seq_len(ng), function(k) {
    fitted(bc, newdata = onec(k), summary = FALSE)[, 1]
  }))
}, 1L)
say("(c) by-level b holds groups 7..12. frmtmb: %s | brms: %s",
    paste(cf, collapse = " "), paste(cb, collapse = " "))
say("(c) agree at %d of 12", sum(cf == cb, na.rm = TRUE))

## (d) (1 | g) + (1 | h), both unseen
fd <- frm(bf(y ~ x + (1 | g) + (1 | h)), family = gaussian(), data = dg)
bd <- bfix(y ~ x + (1 | g) + (1 | h), dg)
ndd <- data.frame(x = 0, g = factor("n1"), h = factor("m1"))
grid <- expand.grid(k = seq_len(ng), l = 1:8)
oned <- function(k, l) data.frame(x = 0, g = factor(k, levels = 1:ng),
                                  h = factor(l, levels = 1:8))
sf <- lapply(seq_len(nrow(grid)), function(i) {
  unname(fitted(fd, newdata = oned(grid$k[i], grid$l[i]))[, 1:2, drop = FALSE])
})
sb <- lapply(seq_len(nrow(grid)), function(i) {
  fitted(bd, newdata = oned(grid$k[i], grid$l[i]), summary = FALSE)[, 1]
})
df_ <- vapply(1:12, function(s) {
  set.seed(s)
  v <- unname(fitted(fd, newdata = ndd, allow_new_levels = TRUE,
                     sample_new_levels = "old_levels")[, 1:2, drop = FALSE])
  i <- which_level(v, sf)
  if (is.na(i)) "NA" else paste0(grid$k[i], ",", grid$l[i])
}, "")
db <- vapply(1:12, function(s) {
  set.seed(s)
  v <- fitted(bd, newdata = ndd, allow_new_levels = TRUE,
              sample_new_levels = "old_levels", summary = FALSE)[, 1]
  i <- which_level(v, sb)
  if (is.na(i)) "NA" else paste0(grid$k[i], ",", grid$l[i])
}, "")
say("(d) (g,h) read. frmtmb: %s", paste(df_, collapse = " "))
say("(d)             brms:   %s", paste(db, collapse = " "))
say("(d) agree at %d of 12", sum(df_ == db))
say("done")
