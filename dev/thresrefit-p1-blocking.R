## Punch round 1. The reviewer's own constructions, on either build.
##   FRMTMB_LIB=base -> 0.64.0 reference build
##   unset           -> the lane library
## Blocking 1 (interior loss, every coding), blocking 2 (influence(data =
## )), finding 3 (hand-written thres(K) on a factor), and the
## influence(groups = ) case that deletes a whole thres(gr = ) level.
arm <- Sys.getenv("FRMTMB_LIB", "lane")
lib <- if (identical(arm, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("ARM", arm, "frmtmb", as.character(packageVersion("frmtmb")), "\n")
say <- function(...) cat(..., "\n", sep = "")
cap <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(list(msg = conditionMessage(e)), class = "probe_err")
    }),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  list(value = v, warnings = w,
       err = if (inherits(v, "probe_err")) v$msg else NULL)
}
show_row <- function(r, i, lab) {
  if (!is.null(r$err)) { say(lab, " -> ERROR: ", r$err); return(invisible()) }
  m <- r$value$fixed
  say(lab, ": dim ", paste(dim(m), collapse = "x"),
      " NA cells ", sum(is.na(m)), " of ", length(m),
      " | all-NA rows ", sum(apply(m, 1, function(z) all(is.na(z)))))
  if (!is.na(i) && i <= nrow(m)) {
    say(lab, ": row ", i, " = ",
        paste(format(m[i, ], digits = 9), collapse = "  "))
  }
  for (ww in r$warnings) say(lab, ": WARNING ", ww)
}

## ==================================================================
## Blocking 1. Seed 901, n = 60, one row in the INTERIOR category 2.
## ==================================================================
mk_interior <- function(seed = 901, n = 60) {
  set.seed(seed)
  x <- rnorm(n)
  cp <- cbind(plogis(0.9 - 0.5 * x), plogis(1.0 - 0.5 * x),
              plogis(2.0 - 0.5 * x))
  y <- 1L + rowSums(runif(n) > cp)
  # exactly ONE row in the interior category 2: the rest are pushed
  # down, and one is planted when the draw produced none
  i2 <- which(y == 2L)
  if (!length(i2)) y[1L] <- 2L else y[i2[-1L]] <- 1L
  data.frame(x = x, y = y)
}
d <- mk_interior()
say("B1: seed 901 table(y) = ", paste(table(factor(d$y, 1:4)), collapse = "/"))
i2 <- which(d$y == 2L)
say("B1: the single interior row is ", i2)

codings <- list(
  integer = function(v) v,
  ordered = function(v) factor(v, levels = 1:4, ordered = TRUE),
  character = function(v) as.character(v))
rows <- list()
for (nm in names(codings)) {
  dd <- data.frame(x = d$x, y = codings[[nm]](d$y))
  fit <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = dd))
  r <- cap(influence(fit, force = TRUE))
  show_row(r, i2, paste0("B1 ", nm))
  if (is.null(r$err)) rows[[nm]] <- r$value$fixed
}
## the reference the pinned refit must reproduce
ref <- suppressWarnings(frm(bf(y | thres(3) ~ x), family = cumulative(),
                            data = d[-i2, , drop = FALSE]))
say("B1 thres(3) reference = ",
    paste(format(frmtmb:::get_coef.frmtmb_fit(ref), digits = 9),
          collapse = "  "))
if (length(rows) > 1L) {
  for (nm in setdiff(names(rows), "integer")) {
    say("B1 ", nm, " identical() to integer, whole table: ",
        identical(rows[[nm]], rows[["integer"]]))
  }
}

## grouped, seed 902: the single interior row sits in level b
mk_gr_interior <- function(seed = 902, n = 96) {
  set.seed(seed)
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.4), b = c(0.9, 1.0, 2.0), c = c(-0.3, 0.9))
  y <- vapply(seq_len(n), function(i) {
    1L + sum(runif(1) > plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
  }, 1L)
  ib <- which(g == "b" & y == 2L)
  if (!length(ib)) {
    y[which(g == "b")[1L]] <- 2L
  } else {
    y[ib[-1L]] <- 1L
  }
  data.frame(x = x, g = g, y = y)
}
dg <- try(mk_gr_interior(), silent = TRUE)
if (!inherits(dg, "try-error")) {
  print(table(dg$g, dg$y))
  ib <- which(dg$g == "b" & dg$y == 2L)
  say("B1g: the single interior row of level b is ", ib)
  grows <- list()
  for (nm in names(codings)) {
    dd <- data.frame(x = dg$x, g = dg$g, y = codings[[nm]](dg$y))
    fit <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                                family = cumulative(), data = dd))
    if (nm == "integer") {
      th <- fit$spec$responses$y$family[["thres"]]
      say("B1g: fitted nthres = ", paste(th[["nthres"]], collapse = "/"),
          " tau_raw length ", length(fit$estimates$tau_raw))
    }
    r <- cap(influence(fit, force = TRUE))
    show_row(r, ib, paste0("B1g ", nm))
    if (is.null(r$err)) grows[[nm]] <- r$value$fixed
  }
  th <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                             family = cumulative(),
                             data = dg))$spec$responses$y$family[["thres"]]
  ds <- dg[-ib, , drop = FALSE]
  ds$k <- as.integer(th[["nthres"]][match(as.character(ds$g),
                                          th[["levels"]])])
  refg <- suppressWarnings(frm(bf(y | thres(k, gr = g) ~ x),
                               family = cumulative(), data = ds))
  say("B1g thres(k, gr = g) reference = ",
      paste(format(frmtmb:::get_coef.frmtmb_fit(refg), digits = 9),
            collapse = "  "))
  for (nm in setdiff(names(grows), "integer")) {
    say("B1g ", nm, " identical() to integer: ",
        identical(grows[[nm]], grows[["integer"]]))
  }
}

## The grouped case above keeps category 2 in levels a and c, so the
## factor keeps all four levels and even the base build is aligned. The
## reviewer's seed 902 needed the category absent GLOBALLY, which is what
## drops the level from the model frame. Constructed here.
mk_gr_global <- function(seed = 904, n = 96) {
  set.seed(seed)
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
  y <- vapply(seq_len(n), function(i) {
    1L + sum(runif(1) > plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
  }, 1L)
  # category 2 nowhere except ONE row, which sits in level a
  y[y == 2L] <- 1L
  y[which(g == "a")[1L]] <- 2L
  data.frame(x = x, g = g, y = y)
}
dh <- mk_gr_global()
print(table(dh$g, dh$y))
ih <- which(dh$y == 2L)
say("B1h: category 2 has ", length(ih), " row(s) in the whole data, row ", ih)
hrows <- list()
for (nm in names(codings)) {
  dd <- data.frame(x = dh$x, g = dh$g, y = codings[[nm]](dh$y))
  fit <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                              family = cumulative(), data = dd))
  if (nm == "integer") {
    th <- fit$spec$responses$y$family[["thres"]]
    say("B1h: fitted nthres = ", paste(th[["nthres"]], collapse = "/"),
        " tau_raw length ", length(fit$estimates$tau_raw))
  }
  r <- cap(influence(fit, force = TRUE))
  show_row(r, ih, paste0("B1h ", nm))
  if (is.null(r$err)) hrows[[nm]] <- r$value$fixed
}
thh <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                            family = cumulative(),
                            data = dh))$spec$responses$y$family[["thres"]]
dsh <- dh[-ih, , drop = FALSE]
dsh$k <- as.integer(thh[["nthres"]][match(as.character(dsh$g),
                                          thh[["levels"]])])
refh <- suppressWarnings(frm(bf(y | thres(k, gr = g) ~ x),
                             family = cumulative(), data = dsh))
say("B1h thres(k, gr = g) reference = ",
    paste(format(frmtmb:::get_coef.frmtmb_fit(refh), digits = 9),
          collapse = "  "))
for (nm in setdiff(names(hrows), "integer")) {
  say("B1h ", nm, " identical() to integer: ",
      identical(hrows[[nm]], hrows[["integer"]]))
}

## ==================================================================
## Finding 3. bf(y | thres(K) ~ x) with an ordered factor.
## ==================================================================
mk_top <- function(seed = 501, n = 50) {
  set.seed(seed)
  x <- rnorm(n)
  cp <- cbind(plogis(-0.7 - 0.5 * x), plogis(0.6 - 0.5 * x),
              plogis(2.3 - 0.5 * x))
  y <- 1L + rowSums(runif(n) > cp)
  top <- which(y == 4L)
  y[top[-1L]] <- 3L
  data.frame(x = x, y = y)
}
dt <- mk_top()
it <- which(dt$y == 4L)
say("F3: seed 501 table = ", paste(table(factor(dt$y, 1:4)), collapse = "/"),
    " top row ", it)
f3 <- list()
for (nm in names(codings)) {
  dd <- data.frame(x = dt$x, y = codings[[nm]](dt$y))
  fit <- suppressWarnings(frm(bf(y | thres(3) ~ x), family = cumulative(),
                              data = dd))
  r <- cap(influence(fit, force = TRUE))
  show_row(r, it, paste0("F3 thres(3) ", nm))
  if (is.null(r$err)) f3[[nm]] <- r$value$fixed
}
for (nm in setdiff(names(f3), "integer")) {
  say("F3 ", nm, " identical() to integer: ",
      identical(f3[[nm]], f3[["integer"]]))
}

## ==================================================================
## Blocking 2. influence(data = ) whose categories differ from the fit.
## ==================================================================
set.seed(1101)
n <- 50
x <- rnorm(n)
y <- 1L + rowSums(cbind(runif(n) > plogis(-0.4 - 0.5 * x),
                        runif(n) > plogis(1.2 - 0.5 * x)))
d1 <- data.frame(x = x, y = pmin(y, 3L))
fit1 <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = d1))
say("B2: fit has n_tau = ", length(fit1$estimates$tau_raw),
    " table = ", paste(table(factor(d1$y, 1:4)), collapse = "/"))
d1b <- d1
d1b$y[1:2] <- 4L
r <- cap(influence(fit1, data = d1b, force = TRUE))
show_row(r, NA, "B2 case 1 (a higher category in data = )")

set.seed(1103)
k <- 96
g <- factor(rep(c("a", "b", "c"), length.out = k))
xx <- rnorm(k)
tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
yy <- vapply(seq_len(k), function(i) {
  1L + sum(runif(1) > plogis(tau[[as.character(g[i])]] - 0.5 * xx[i]))
}, 1L)
dg2 <- data.frame(x = xx, g = g, y = yy)
fitg <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                             family = cumulative(), data = dg2))
say("B2: grouped fit nthres = ",
    paste(fitg$spec$responses$y$family$thres$nthres, collapse = "/"))
dg2b <- dg2
levels(dg2b$g) <- c(levels(dg2$g), "d")
dg2b$g[1:4] <- "d"
r <- cap(influence(fitg, data = dg2b, force = TRUE))
show_row(r, NA, "B2 case 2 (a new thres(gr = ) level in data = )")

## and influence(groups = ) that deletes a whole thres(gr = ) level
dg3 <- dg2
dg3$idg <- dg3$g
fitg3 <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x + (1 | idg)),
                              family = cumulative(), data = dg3))
r <- cap(influence(fitg3, groups = "idg"))
show_row(r, NA, "B2 case 3 (groups = deletes a whole thres(gr = ) level)")

## the control: a partly failed table must warn, not stop, and a table
## with nothing wrong must stay silent
r <- cap(influence(fit1, force = TRUE))
say("B2 control (nothing wrong): warnings = ", length(r$warnings),
    " NA cells = ", if (is.null(r$err)) sum(is.na(r$value$fixed)) else NA)
