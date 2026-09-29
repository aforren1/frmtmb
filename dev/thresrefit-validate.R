## Validation for the threshold-pinning change. Runs on either build:
##   FRMTMB_LIB=base  -> the 0.64.0 reference build (the "before" arm)
##   unset            -> the lane library (the "after" arm)
arm <- Sys.getenv("FRMTMB_LIB", "lane")
lib <- if (identical(arm, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("ARM", arm, "frmtmb", as.character(packageVersion("frmtmb")), "\n")
say <- function(...) cat(..., "\n", sep = "")
try_msg <- function(expr) {
  tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
}
relerr <- function(a, b) abs(a - b) / max(abs(a), abs(b))

## ==================================================================
## 1. influence() on an ordinal fit whose top category has ONE row
## ==================================================================
mk1 <- function(seed = 501, n = 50) {
  set.seed(seed)
  x <- rnorm(n)
  cp <- cbind(plogis(-0.7 - 0.5 * x), plogis(0.6 - 0.5 * x),
              plogis(2.3 - 0.5 * x))
  y <- 1L + rowSums(runif(n) > cp)
  ## keep exactly one row in the top category, so one deletion removes it
  top <- which(y == 4L)
  if (length(top) > 1L) y[top[-1L]] <- 3L
  if (!length(top)) stop("seed reached no top category")
  data.frame(x = x, y = y)
}
d1 <- mk1()
say("1: seed 501 n = 50 table(y) = ",
    paste(table(factor(d1$y, 1:4)), collapse = "/"))
itop <- which(d1$y == 4L)
say("1: the single top-category row is ", itop)
for (famnm in c("cumulative", "sratio")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  fit <- suppressWarnings(frm(bf(y ~ x), family = f, data = d1))
  inf <- try_msg(suppressWarnings(influence(fit, force = TRUE)))
  if (is.character(inf)) { say("1 ", famnm, ": influence -> ", inf); next }
  say("1 ", famnm, ": colnames = ",
      paste(colnames(inf$fixed), collapse = ","))
  say("1 ", famnm, ": NA cells in the whole table = ", sum(is.na(inf$fixed)))
  say("1 ", famnm, ": row ", itop, " = ",
      paste(format(inf$fixed[itop, ], digits = 10), collapse = " "))
  cd <- cooks.distance(inf)
  say("1 ", famnm, ": cooks.distance NA count = ", sum(is.na(cd)),
      "; value at row ", itop, " = ", format(cd[itop], digits = 8))
  ## the reference the pinned refit must reproduce: thres(3) on the
  ## deleted data, fit from scratch
  dsub <- d1[-itop, , drop = FALSE]
  fref <- suppressWarnings(frm(bf(y | thres(3) ~ x), family = f,
                               data = dsub))
  cref <- c(coef_x = as.numeric(fref$estimates$beta),
            fref$estimates$tau_raw)
  got <- as.numeric(inf$fixed[itop, ])
  say("1 ", famnm, ": thres(3) reference = ",
      paste(format(cref, digits = 10), collapse = " "))
  if (!anyNA(got)) {
    say("1 ", famnm, ": max relative difference to the reference = ",
        format(max(mapply(relerr, got, cref)), digits = 4))
  }
}

## ==================================================================
## 2. influence() with GROUPED thresholds
## ==================================================================
mk2 <- function(seed = 502, n = 96) {
  set.seed(seed)
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
  y <- integer(n)
  for (i in seq_len(n)) {
    tt <- tau[[as.character(g[i])]]
    y[i] <- 1L + sum(runif(1) > plogis(tt - 0.5 * x[i]))
  }
  ## level a keeps exactly one row in its top category
  ia <- which(g == "a" & y == 4L)
  if (length(ia) > 1L) y[ia[-1L]] <- 3L
  if (!length(ia)) stop("seed reached no top category in level a")
  data.frame(x = x, g = g, y = y)
}
d2 <- try_msg(mk2())
if (is.character(d2)) say("2: ", d2) else {
  print(table(d2$g, d2$y))
  fg <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                            family = cumulative(), data = d2))
  say("2: fitted tau_raw length = ", length(fg$estimates$tau_raw),
      " per-level counts = ",
      paste(fg$spec$responses$y$family$thres$nthres, collapse = ","))
  ia <- which(d2$g == "a" & d2$y == 4L)
  say("2: the single top row of level a is ", ia)
  infg <- try_msg(suppressWarnings(influence(fg, force = TRUE)))
  if (is.character(infg)) say("2: influence -> ", infg) else {
    say("2: colnames = ", paste(colnames(infg$fixed), collapse = ","))
    say("2: NA cells = ", sum(is.na(infg$fixed)))
    say("2: row ", ia, " = ",
        paste(format(infg$fixed[ia, ], digits = 8), collapse = " "))
    cdg <- cooks.distance(infg)
    say("2: cooks.distance NA count = ", sum(is.na(cdg)),
        "; value at row ", ia, " = ", format(cdg[ia], digits = 8))
    ## reference: thres(k, gr = g) with the FITTED counts on the
    ## deleted data
    dsub <- d2[-ia, , drop = FALSE]
    np <- fg$spec$responses$y$family$thres$nthres
    lv <- fg$spec$responses$y$family$thres$levels
    dsub$k <- as.integer(np[match(as.character(dsub$g), lv)])
    fref <- suppressWarnings(frm(bf(y | thres(k, gr = g) ~ x),
                                family = cumulative(), data = dsub))
    cref <- c(as.numeric(fref$estimates$beta), fref$estimates$tau_raw)
    got <- as.numeric(infg$fixed[ia, ])
    say("2: reference = ", paste(format(cref, digits = 8), collapse = " "))
    if (!anyNA(got) && length(got) == length(cref)) {
      say("2: max relative difference = ",
          format(max(mapply(relerr, got, cref)), digits = 4))
    }
  }
}

## ==================================================================
## 3. frm_bootstrap(): the layout must be constant, and the replicate
##    that lost the category must match a thres(K) pinned fit
## ==================================================================
mk3 <- function(seed = 202, n = 40) {
  set.seed(seed)
  x <- rnorm(n)
  cp <- cbind(plogis(-0.6 - 0.5 * x), plogis(0.5 - 0.5 * x),
              plogis(2.6 - 0.5 * x))
  data.frame(x = x, y = 1L + rowSums(runif(n) > cp))
}
d3 <- mk3()
say("3: seed 202 n = 40 table(y) = ",
    paste(table(factor(d3$y, 1:4)), collapse = "/"))
FUNt <- function(ff) c(fixef(ff, flatten = TRUE), tau = ff$estimates$tau_raw)
for (famnm in c("cumulative", "sratio")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  fit <- suppressWarnings(frm(bf(y ~ x), family = f, data = d3))
  set.seed(11)
  sims <- simulate(fit, nsim = 60, re_formula = NA)
  lost <- which(vapply(sims, function(v) max(as.integer(v)), 1L) < 4L)
  say("3 ", famnm, ": replicates without the top category (seed 11, 60) = ",
      paste(lost, collapse = ","))
  bs <- try_msg(suppressWarnings(frm_bootstrap(fit, FUN = FUNt, nsim = 60,
                                               seed = 11)))
  if (is.character(bs)) { say("3 ", famnm, ": bootstrap -> ", bs); next }
  say("3 ", famnm, ": t dim = ", paste(dim(bs$t), collapse = "x"),
      " colnames = ", paste(colnames(bs$t), collapse = ","),
      " NA cells = ", sum(is.na(bs$t)))
  for (b in lost) {
    yb <- as.integer(sims[[b]])
    db <- data.frame(x = d3$x, y = yb)
    fp <- suppressWarnings(frm(bf(y | thres(3) ~ x), family = f, data = db))
    rf <- suppressWarnings(refit(fit, yb))
    say("3 ", famnm, ": replicate ", b,
        " refit logLik = ", format(as.numeric(logLik(rf)), digits = 15),
        " thres(3) logLik = ", format(as.numeric(logLik(fp)), digits = 15),
        " relative difference = ",
        format(relerr(as.numeric(logLik(rf)), as.numeric(logLik(fp))),
               digits = 4))
    say("3 ", famnm, ": replicate ", b, " tau_raw length refit = ",
        length(rf$estimates$tau_raw), " pinned = ",
        length(fp$estimates$tau_raw))
  }
}

## grouped bootstrap
d2b <- try_msg(mk2(503))
if (!is.character(d2b)) {
  fg <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                            family = cumulative(), data = d2b))
  np <- fg$spec$responses$y$family$thres$nthres
  lv <- fg$spec$responses$y$family$thres$levels
  say("3g: fitted per-level counts = ", paste(np, collapse = ","))
  set.seed(13)
  sg <- simulate(fg, nsim = 40, re_formula = NA)
  lostg <- which(vapply(sg, function(v) {
    any(tapply(as.integer(v), d2b$g, max) - 1L < np)
  }, TRUE))
  say("3g: replicates where a level lost its top category = ",
      length(lostg), " (", paste(utils::head(lostg, 8), collapse = ","), ")")
  bsg <- try_msg(suppressWarnings(
    frm_bootstrap(fg, FUN = FUNt, nsim = 40, seed = 13)))
  if (is.character(bsg)) say("3g: bootstrap -> ", bsg) else {
    say("3g: t dim = ", paste(dim(bsg$t), collapse = "x"),
        " NA cells = ", sum(is.na(bsg$t)))
  }
  if (length(lostg)) {
    b <- lostg[1L]
    yb <- as.integer(sg[[b]])
    dbb <- data.frame(x = d2b$x, g = d2b$g, y = yb,
                      k = as.integer(np[match(as.character(d2b$g), lv)]))
    fp <- suppressWarnings(frm(bf(y | thres(k, gr = g) ~ x),
                              family = cumulative(), data = dbb))
    rf <- suppressWarnings(refit(fg, yb))
    say("3g: replicate ", b, " refit logLik = ",
        format(as.numeric(logLik(rf)), digits = 15),
        " pinned logLik = ", format(as.numeric(logLik(fp)), digits = 15),
        " relative difference = ",
        format(relerr(as.numeric(logLik(rf)), as.numeric(logLik(fp))),
               digits = 4))
  }
}

## ==================================================================
## 4. categorical(): the same shape of defect, or not
## ==================================================================
mk4 <- function(seed = 504, n = 50) {
  set.seed(seed)
  x <- rnorm(n)
  P <- cbind(1, exp(0.2 + 0.4 * x), exp(-2.6 + 0.3 * x))
  P <- P / rowSums(P)
  y <- apply(P, 1, function(p) sample.int(3L, 1L, prob = p))
  top <- which(y == 3L)
  if (length(top) > 1L) y[top[-1L]] <- 2L
  data.frame(x = x, y = factor(y))
}
d4 <- mk4()
say("4: table(y) = ", paste(table(d4$y), collapse = "/"))
f4 <- suppressWarnings(frm(bf(y ~ x), family = categorical(), data = d4))
say("4: dpars = ", paste(names(f4$spec$responses$y$dpars), collapse = ","))
i3 <- which(d4$y == "3")
say("4: the single level-3 row is ", i3)
i4 <- try_msg(suppressWarnings(influence(f4, force = TRUE)))
if (is.character(i4)) say("4: influence -> ", i4) else {
  say("4: colnames = ", paste(colnames(i4$fixed), collapse = ","))
  say("4: NA cells = ", sum(is.na(i4$fixed)))
  say("4: row ", i3, " = ",
      paste(format(i4$fixed[i3, ], digits = 6), collapse = " "))
}
bs4 <- try_msg(suppressWarnings(frm_bootstrap(f4, nsim = 20, seed = 21)))
if (is.character(bs4)) say("4: bootstrap -> ", bs4) else {
  say("4: boot t dim = ", paste(dim(bs4$t), collapse = "x"),
      " NA cells = ", sum(is.na(bs4$t)))
}

## ==================================================================
## 5. B: the prior draw on a threshold vector, both directions
## ==================================================================
## (a) more than one threshold: no draw may recycle
db <- data.frame(x = rnorm(40), y = rep(1:4, 10))
for (famnm in c("cratio", "acat", "sratio", "cumulative")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  r <- try_msg(frm_simulate(bf(y ~ x), family = f, data = db,
                            prior = set_prior("normal(0, 2)",
                                              class = "Intercept") +
                              set_prior("normal(0, 1)", class = "b"),
                            nsim = 2, seed = 5))
  say("5a ", famnm, " (3 thresholds) -> ",
      if (is.character(r)) r else "DREW (see pars below)")
  if (!is.character(r)) print(attr(r, "pars"))
}
## (b) exactly one threshold: the draw is one number for one parameter,
##     so there is nothing to recycle and it must still work
db2 <- data.frame(x = rnorm(40), y = rep(1:2, 20))
for (famnm in c("cratio", "acat", "sratio", "cumulative")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  r <- try_msg(frm_simulate(bf(y ~ x), family = f, data = db2,
                            prior = set_prior("normal(0, 2)",
                                              class = "Intercept") +
                              set_prior("normal(0, 1)", class = "b"),
                            nsim = 2, seed = 5))
  if (is.character(r)) say("5b ", famnm, " (1 threshold) -> ", r) else {
    p <- attr(r, "pars")
    say("5b ", famnm, " (1 threshold) -> drew ",
        paste(names(p), collapse = ","), " = ",
        paste(format(unlist(p[1, ]), digits = 6), collapse = " "))
  }
}
## (c) the internal draw itself, on an unordered 3-threshold entry
e <- list(comp = "tau_raw", idx = 1:3, scale = "internal",
          dist = frmtmb:::prior_normal(0, 2))
r <- try_msg(frmtmb:::draw_prior_entry(e, "tau_raw[1:3]"))
say("5c draw_prior_entry on idx 1:3 -> ",
    if (is.character(r)) r else paste("returned", format(r, digits = 8)))
est <- list(tau_raw = c(0, 0, 0))
r2 <- try_msg(frmtmb:::draw_prior_pars(est, list(e), "tau_raw[1:3]"))
say("5c draw_prior_pars -> ",
    if (is.character(r2)) r2 else
      paste("tau_raw =", paste(format(r2$est$tau_raw, digits = 8),
                               collapse = " ")))
