## The agreement between a pinned leave-one-out refit and a fit that
## pins the count by hand, reported per coefficient rather than as one
## max that the unidentified threshold dominates.
.libPaths(c("C:/Users/adf44/source/r/wt-thresrefit-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
say <- function(...) cat(..., "\n", sep = "")

one_top <- function(seed = 501, n = 50) {
  set.seed(seed)
  x <- rnorm(n)
  cp <- cbind(plogis(-0.7 - 0.5 * x), plogis(0.6 - 0.5 * x),
              plogis(2.3 - 0.5 * x))
  y <- 1L + rowSums(runif(n) > cp)
  top <- which(y == 4L)
  y[top[-1L]] <- 3L
  data.frame(x = x, y = y)
}
dd <- one_top()
itop <- which(dd$y == 4L)
for (famnm in c("cumulative", "sratio")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  fit <- suppressWarnings(frm(bf(y ~ x), family = f, data = dd))
  inf <- suppressWarnings(influence(fit, force = TRUE))
  ref <- suppressWarnings(frm(bf(y | thres(3) ~ x), family = f,
                              data = dd[-itop, , drop = FALSE]))
  want <- frmtmb:::get_coef.frmtmb_fit(ref)
  got <- inf$fixed[itop, ]
  rel <- abs(got - want) / abs(want)
  say("== ", famnm)
  for (nm in names(want)) {
    say("  ", nm, ": refit ", format(got[[nm]], digits = 15),
        " pinned ", format(want[[nm]], digits = 15),
        " relative ", format(rel[[nm]], digits = 4))
  }
  sp <- apply(inf$fixed, 2, stats::sd)
  say("  |difference| / sd(column): ",
      paste(paste0(names(want), "=",
                   format(abs(got - want) / sp, digits = 3)),
            collapse = " "))
}

## grouped
grp <- function(seed = 502, n = 96) {
  set.seed(seed)
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
  y <- vapply(seq_len(n), function(i) {
    1L + sum(runif(1) > plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
  }, 1L)
  ia <- which(g == "a" & y == 4L)
  y[ia[-1L]] <- 3L
  data.frame(x = x, g = g, y = y)
}
dg <- grp()
ia <- which(dg$g == "a" & dg$y == 4L)
fg <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                          family = cumulative(), data = dg))
th <- fg$spec$responses$y$family[["thres"]]
inf <- suppressWarnings(influence(fg, force = TRUE))
ds <- dg[-ia, , drop = FALSE]
ds$k <- as.integer(th[["nthres"]][match(as.character(ds$g), th[["levels"]])])
ref <- suppressWarnings(frm(bf(y | thres(k, gr = g) ~ x),
                           family = cumulative(), data = ds))
want <- frmtmb:::get_coef.frmtmb_fit(ref)
got <- inf$fixed[ia, ]
say("== grouped cumulative, per-level counts ",
    paste(th[["nthres"]], collapse = "/"))
for (nm in names(want)) {
  say("  ", nm, ": refit ", format(got[[nm]], digits = 15),
      " pinned ", format(want[[nm]], digits = 15),
      " relative ", format(abs(got[[nm]] - want[[nm]]) / abs(want[[nm]]),
                           digits = 4))
}
sp <- apply(inf$fixed, 2, stats::sd)
say("  |difference| / sd(column): ",
    paste(paste0(names(want), "=", format(abs(got - want) / sp, digits = 3)),
          collapse = " "))
