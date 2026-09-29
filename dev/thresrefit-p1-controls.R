## Punch round 1, the ABSENT cases for the two new refusals: what must
## NOT be refused, and the partial-failure warning rather than the stop.
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

## ------------------------------------------------------------------
## Control 1. influence(groups = ) on a thres(gr = ) model where the
## deleted factor is NOT the threshold grouping factor. Every level of
## thres(gr = ) keeps rows, so nothing may be refused.
## ------------------------------------------------------------------
set.seed(1201)
n <- 120
g <- factor(rep(c("a", "b", "c"), length.out = n))
id <- factor(rep(1:8, length.out = n))
x <- rnorm(n)
tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
y <- vapply(seq_len(n), function(i) {
  1L + sum(runif(1) > plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
}, 1L)
dd <- data.frame(x = x, g = g, id = id, y = y)
print(table(dd$g, dd$y))
fit <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x + (1 | id)),
                            family = cumulative(), data = dd))
say("C1: fitted nthres = ",
    paste(fit$spec$responses$y$family$thres$nthres, collapse = "/"),
    " tau_raw length ", length(fit$estimates$tau_raw))
say("C1: rows of g surviving each id deletion (min over id and level) = ",
    min(table(dd$g, dd$id) %*% (1 - diag(8)) ))
r <- cap(influence(fit, groups = "id"))
if (!is.null(r$err)) say("C1 -> ERROR: ", r$err) else {
  m <- r$value$fixed
  say("C1: dim ", paste(dim(m), collapse = "x"), " NA cells ",
      sum(is.na(m)), " of ", length(m), " warnings ", length(r$warnings))
  say("C1: cooks.distance NA count = ",
      sum(is.na(cooks.distance(r$value))))
}

## ------------------------------------------------------------------
## Control 2. A PARTIAL failure must warn and return the table, not
## stop: exactly one row of `data` reaches a category the fit never saw,
## so the single deletion that removes it succeeds and the rest fail.
## ------------------------------------------------------------------
set.seed(1202)
n2 <- 40
x2 <- rnorm(n2)
y2 <- 1L + rowSums(cbind(runif(n2) > plogis(-0.4 - 0.5 * x2),
                         runif(n2) > plogis(1.2 - 0.5 * x2)))
d2 <- data.frame(x = x2, y = pmin(y2, 3L))
fit2 <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = d2))
say("C2: fit n_tau = ", length(fit2$estimates$tau_raw),
    " table = ", paste(table(factor(d2$y, 1:4)), collapse = "/"))
d2b <- d2
d2b$y[7L] <- 4L
r <- cap(influence(fit2, data = d2b, force = TRUE))
if (!is.null(r$err)) say("C2 -> ERROR: ", r$err) else {
  m <- r$value$fixed
  say("C2: dim ", paste(dim(m), collapse = "x"), " NA cells ",
      sum(is.na(m)), " of ", length(m),
      " | all-NA rows ", sum(apply(m, 1, function(z) all(is.na(z)))))
  say("C2: the ONE row that refits is ",
      paste(which(apply(m, 1, function(z) !any(is.na(z)))), collapse = ","))
  for (w in r$warnings) say("C2: WARNING ", w)
}

## ------------------------------------------------------------------
## Control 3. The pin must stay inert: a plain ordinal fit with nothing
## to lose, and four non-ordinal families, must warn about nothing.
## ------------------------------------------------------------------
set.seed(1203)
n3 <- 40
d3 <- data.frame(x = rnorm(n3))
d3$y <- 1L + rowSums(cbind(runif(n3) > plogis(-0.8 - 0.5 * d3$x),
                           runif(n3) > plogis(0.8 - 0.5 * d3$x)))
say("C3 ordinal table = ", paste(table(factor(d3$y, 1:3)), collapse = "/"))
for (nm in c("cumulative", "sratio", "cratio", "acat")) {
  f <- get(nm, envir = asNamespace("frmtmb"))()
  fit3 <- suppressWarnings(frm(bf(y ~ x), family = f, data = d3))
  r <- cap(influence(fit3, force = TRUE))
  say("C3 ", nm, ": NA cells ",
      if (is.null(r$err)) sum(is.na(r$value$fixed)) else NA,
      " warnings ", length(r$warnings),
      if (!is.null(r$err)) paste0(" ERROR: ", r$err) else "")
}
d4 <- data.frame(x = rnorm(40))
d4$yg <- rnorm(40, d4$x)
d4$yb <- rbinom(40, 1, plogis(d4$x))
d4$yp <- rpois(40, exp(0.5 + 0.3 * d4$x))
for (nm in c("gaussian", "bernoulli", "poisson")) {
  f <- get(nm, envir = asNamespace("frmtmb"))()
  yv <- switch(nm, gaussian = "yg", bernoulli = "yb", poisson = "yp")
  fm <- stats::as.formula(paste(yv, "~ x"))
  fit4 <- suppressWarnings(frm(bf(fm), family = f, data = d4))
  r <- cap(influence(fit4, force = TRUE))
  say("C3 ", nm, ": NA cells ",
      if (is.null(r$err)) sum(is.na(r$value$fixed)) else NA,
      " warnings ", length(r$warnings),
      " pin ", if (is.null(frmtmb:::thres_pin_of_fit(fit4))) "NULL" else "set")
}

## ------------------------------------------------------------------
## Control 4. The accepts_aterms half of the gate, which no family in
## the package makes false: an ordinal family that does not declare
## thres must not be pinned.
## ------------------------------------------------------------------
fam_no_thres <- frmtmb_family("ordy", type = "ordinal",
                              accepts_aterms = "weights",
                              dpars = "mu", links = list(mu = "identity"),
                              lpdf = function(y, dpars, aterms) 0 * y)
shim <- list(frame = list(par_template = list(tau_raw = c(0, 0, 0)),
                          y_levels = list(y = c("1", "2", "3", "4"))),
             spec = list(responses = list(
               y = list(resp_name = "y", family = fam_no_thres))))
say("C4: ordinal family without `thres` in accepts_aterms -> pin ",
    if (is.null(frmtmb:::thres_pin_of_fit(shim))) "NULL" else "SET")
fam_thres <- frmtmb_family("ordy2", type = "ordinal",
                           accepts_aterms = c("weights", "thres"),
                           dpars = "mu", links = list(mu = "identity"),
                           lpdf = function(y, dpars, aterms) 0 * y)
shim2 <- shim
shim2$spec$responses$y$family <- fam_thres
p <- frmtmb:::thres_pin_of_fit(shim2)
say("C4: the same family WITH `thres` -> pin ",
    if (is.null(p)) "NULL" else paste0("nthres=", p$y$nthres,
                                       " levels_y=",
                                       paste(p$y$levels_y, collapse = ",")))
