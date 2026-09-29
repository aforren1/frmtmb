## RE-CHECK item 2: a factor response whose category labels are NOT
## "1".."K", and an ordered factor with a level that no row takes even in
## the FULL fit. The recode maps codes through labels, so a label set that
## is not sorted, not numeric, or wider than the data is where it can
## misread. Checked on the FULL fit and on the deletion refits.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, "\n\n")
options(width = 160, digits = 17)

## integer 1..4, exactly one row in the interior category 2 and one in the
## top category 4, so both kinds of emptying are exercised at once
mk <- function(seed = 901, n = 60) {
  set.seed(seed)
  x <- stats::rnorm(n)
  cp <- cbind(stats::plogis(-0.7 - 0.5 * x), stats::plogis(0.6 - 0.5 * x),
              stats::plogis(2.0 - 0.5 * x))
  y <- 1L + rowSums(stats::runif(n) > cp)
  mid <- which(y == 2L)
  y[mid[-1L]] <- 1L
  top <- which(y == 4L)
  y[top[-1L]] <- 3L
  data.frame(x = x, y = y)
}
dd <- mk()
cat("table(y) =", paste(table(dd$y), collapse = "/"),
    "  interior row =", which(dd$y == 2L),
    "  top row =", which(dd$y == 4L), "\n\n")

LAB <- c("lo", "mid", "hi", "top")
codings <- list(
  integer   = function(y) y,
  ordered   = function(y) factor(y, levels = 1:4, ordered = TRUE),
  ## labels that are not numbers, declared in the intended order; the
  ## alphabetical order would be hi < lo < mid < top, so a recode that
  ## used sort() instead of the fitted labels would scramble this
  labelled  = function(y) factor(LAB[y], levels = LAB, ordered = TRUE),
  ## labels declared in DESCENDING alphabetical order, so sorting is wrong
  descend   = function(y) factor(c("d", "c", "b", "a")[y],
                                levels = c("d", "c", "b", "a"),
                                ordered = TRUE),
  ## five declared levels, the fifth taken by no row even in the FULL fit
  wide      = function(y) factor(y, levels = 1:5, ordered = TRUE),
  character = function(y) as.character(y),
  ## a character vector of non-numeric labels: what does the family do?
  charlab   = function(y) LAB[y]
)

res <- list()
for (nm in names(codings)) {
  d <- dd
  d$y <- codings[[nm]](dd$y)
  fit <- try(suppressWarnings(frm(bf(y ~ x), family = cumulative(),
                                  data = d)), silent = TRUE)
  if (inherits(fit, "try-error")) {
    cat("==", nm, ": frm() REFUSED:",
        substr(gsub("\n", " ", conditionMessage(attr(fit, "condition"))),
               1, 120), "\n\n")
    next
  }
  lvf <- fit$frame[["y_levels"]][["y"]]
  cat("==", nm, ": n_tau =", length(fit$estimates[["tau_raw"]]),
      " y_levels =", if (is.null(lvf)) "NULL" else
        paste(lvf, collapse = ","), "\n")
  cat("   FULL fit logLik =", sprintf("%.14f", logLik(fit)),
      " coefs =", paste(sprintf("%.14g",
                                frmtmb:::get_coef.frmtmb_fit(fit)),
                        collapse = " "), "\n")
  if (exists("thres_pin_of_fit", envir = asNamespace("frmtmb"),
             inherits = FALSE)) {
    p <- frmtmb:::thres_pin_of_fit(fit)
    cat("   pin: nthres =", p$y$nthres, " levels_y =",
        if (is.null(p$y$levels_y)) "NULL" else
          paste(p$y$levels_y, collapse = ","), "\n")
  }
  nw <- 0L
  inf <- withCallingHandlers(try(influence(fit, force = TRUE), silent = TRUE),
                             warning = function(w) {
                               nw <<- nw + 1L
                               cat("     WARNING:",
                                   gsub("\n", " ", conditionMessage(w)), "\n")
                               invokeRestart("muffleWarning")
                             })
  if (inherits(inf, "try-error")) {
    cat("   influence() ERROR:",
        gsub("\n", " ", conditionMessage(attr(inf, "condition"))), "\n\n")
    next
  }
  res[[nm]] <- inf$fixed
  cat("   influence: NA cells =", sum(is.na(inf$fixed)),
      " warnings =", nw,
      " NA in interior row =", sum(is.na(inf$fixed[which(dd$y == 2L), ])),
      " NA in top row =", sum(is.na(inf$fixed[which(dd$y == 4L), ])), "\n")
  cat("   interior row:",
      paste(signif(inf$fixed[which(dd$y == 2L), ], 9), collapse = "  "), "\n")
  cat("   top row     :",
      paste(signif(inf$fixed[which(dd$y == 4L), ], 9), collapse = "  "), "\n")
  cat("\n")
}
cat("#### identical() against the integer coding, whole table\n")
if ("integer" %in% names(res)) {
  for (nm in setdiff(names(res), "integer")) {
    cat("  ", nm, ":", identical(res[[nm]], res[["integer"]]),
        "  dim match:", identical(dim(res[[nm]]), dim(res[["integer"]])),
        "\n")
  }
}

## the hand-written thres(K) variant on the labelled coding
cat("\n#### bf(y | thres(3) ~ x) on the labelled ordered factor\n")
d <- dd
d$y <- codings$labelled(dd$y)
f2 <- try(suppressWarnings(frm(bf(y | thres(3) ~ x), family = cumulative(),
                               data = d)), silent = TRUE)
if (inherits(f2, "try-error")) {
  cat("   frm() REFUSED:",
      gsub("\n", " ", conditionMessage(attr(f2, "condition"))), "\n")
} else {
  nw <- 0L
  i2 <- withCallingHandlers(try(influence(f2, force = TRUE), silent = TRUE),
                            warning = function(w) {
                              nw <<- nw + 1L; invokeRestart("muffleWarning")
                            })
  if (inherits(i2, "try-error")) {
    cat("   influence() ERROR:",
        gsub("\n", " ", conditionMessage(attr(i2, "condition"))), "\n")
  } else {
    cat("   NA cells =", sum(is.na(i2$fixed)), " warnings =", nw,
        "  identical to the integer thres(3) table:",
        {
          di <- dd
          fi <- suppressWarnings(frm(bf(y | thres(3) ~ x),
                                     family = cumulative(), data = di))
          ii <- suppressWarnings(influence(fi, force = TRUE))
          identical(i2$fixed, ii$fixed)
        }, "\n")
  }
}

## grouped, with labelled response levels
cat("\n#### grouped thres(gr = g) with labelled response levels\n")
set.seed(904)
n <- 96
g <- factor(rep(c("a", "b", "c"), length.out = n))
x <- stats::rnorm(n)
tau <- list(a = c(-0.6, 0.6, 2.0), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
yg <- vapply(seq_len(n), function(i) {
  1L + sum(stats::runif(1) >
             stats::plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
}, 1L)
mid <- which(yg == 2L)
yg[mid[-1L]] <- 1L
dg <- data.frame(x = x, g = g, y = yg)
out <- list()
for (nm in c("integer", "labelled")) {
  d <- dg
  d$y <- codings[[nm]](dg$y)
  fit <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                             family = cumulative(), data = d))
  nw <- 0L
  inf <- withCallingHandlers(suppressWarnings(influence(fit, force = TRUE)),
                             warning = function(w) {
                               nw <<- nw + 1L; invokeRestart("muffleWarning")
                             })
  out[[nm]] <- inf$fixed
  cat("  ", nm, ": nthres =",
      paste(fit$spec$responses$y$family[["thres"]][["nthres"]],
            collapse = "/"),
      " NA =", sum(is.na(inf$fixed)), " warnings =", nw, "\n")
  cat("       interior row:",
      paste(signif(inf$fixed[which(dg$y == 2L)[1L], ], 9), collapse = "  "),
      "\n")
}
if (length(out) == 2L) {
  cat("   identical(labelled, integer):",
      identical(out$labelled, out$integer), "\n")
}
cat("DONE rev-21\n")
