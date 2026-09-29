## RE-CHECK item 1: the reviewer's own interior constructions, across
## integer, ordered-factor and character codings, against the
## thres(3) / thres(k, gr = ) reference on the same subset.
## Seeds 901 and 902 are the reviewer's from rev-02; 904 is the global
## interior loss the worker says is needed to empty a level.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, " ver =",
    as.character(utils::packageVersion("frmtmb")), "\n\n")
options(width = 160, digits = 17)

recode <- function(y, how) switch(how,
  integer = y,
  ordered = factor(y, levels = sort(unique(y)), ordered = TRUE),
  character = as.character(y))

## reviewer's rev-02 generator, verbatim
mk <- function(seed = 901, n = 60) {
  set.seed(seed)
  x <- stats::rnorm(n)
  cp <- cbind(stats::plogis(-0.7 - 0.5 * x), stats::plogis(0.6 - 0.5 * x),
              stats::plogis(2.0 - 0.5 * x))
  y <- 1L + rowSums(stats::runif(n) > cp)
  mid <- which(y == 2L)
  y[mid[-1L]] <- 1L
  data.frame(x = x, y = y)
}
mkg <- function(seed = 902, n = 96) {
  set.seed(seed)
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- stats::rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.0), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
  y <- vapply(seq_len(n), function(i) {
    1L + sum(stats::runif(1) >
               stats::plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
  }, 1L)
  mid <- which(y == 2L)
  y[mid[-1L]] <- 1L
  data.frame(x = x, g = g, y = y)
}

arm <- function(dd, irow, grouped, tag, fitted_nthres = NULL) {
  cat("####", tag, "  deleted row =", irow, "\n")
  cat("   table(y) =", paste(table(dd$y), collapse = "/"),
      if (grouped) paste0("  per-level max = ",
                          paste(tapply(dd$y, dd$g, max), collapse = "/")),
      "\n")
  form <- if (grouped) bf(y | thres(gr = g) ~ x) else bf(y ~ x)
  keep <- list()
  warn <- list()
  for (how in c("integer", "ordered", "character")) {
    d <- dd
    d$y <- recode(dd$y, how)
    fit <- try(suppressWarnings(frm(form, family = cumulative(), data = d)),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      cat("  ", how, ": frm() REFUSED\n"); next
    }
    nw <- 0L
    inf <- withCallingHandlers(
      try(influence(fit, force = TRUE), silent = TRUE),
      warning = function(w) {
        nw <<- nw + 1L
        cat("      WARNING:", gsub("\n", " ", conditionMessage(w)), "\n")
        invokeRestart("muffleWarning")
      })
    if (inherits(inf, "try-error")) {
      cat("  ", how, ": influence() ERROR:",
          gsub("\n", " ", conditionMessage(attr(inf, "condition"))), "\n")
      next
    }
    keep[[how]] <- inf$fixed
    warn[[how]] <- nw
    cat("  ", how, ": n_tau =", length(fit$estimates[["tau_raw"]]),
        " ncol =", ncol(inf$fixed), " NA cells =", sum(is.na(inf$fixed)),
        " NA in deleted row =", sum(is.na(inf$fixed[irow, ])),
        " warnings =", nw, "\n")
    cat("       row:", paste(signif(inf$fixed[irow, ], 9), collapse = "  "),
        "\n")
    cat("       cooks NA =", sum(is.na(suppressWarnings(
      cooks.distance(inf)))),
        "  cooks[row] =", suppressWarnings(cooks.distance(inf))[irow], "\n")
  }
  if ("integer" %in% names(keep)) {
    for (how in setdiff(names(keep), "integer")) {
      cat("   identical(", how, ", integer) over the WHOLE table:",
          identical(keep[[how]], keep[["integer"]]), "\n")
    }
  }
  ## the reference: the same subset, integer codes, count pinned by hand
  sub <- dd[-irow, , drop = FALSE]
  if (grouped) {
    sub$k <- as.integer(fitted_nthres[match(as.character(sub$g),
                                            names(fitted_nthres))])
    ref <- try(suppressWarnings(frm(bf(y | thres(k, gr = g) ~ x),
                                    family = cumulative(), data = sub)),
               silent = TRUE)
  } else {
    ref <- try(suppressWarnings(frm(bf(y | thres(3) ~ x),
                                    family = cumulative(), data = sub)),
               silent = TRUE)
  }
  if (inherits(ref, "try-error")) {
    cat("   reference REFUSED:",
        gsub("\n", " ", conditionMessage(attr(ref, "condition"))), "\n\n")
    return(invisible(NULL))
  }
  w <- frmtmb:::get_coef.frmtmb_fit(ref)
  cat("   reference:", paste(signif(w, 9), collapse = "  "), "\n")
  if ("integer" %in% names(keep)) {
    got <- keep[["integer"]][irow, ]
    sp <- apply(keep[["integer"]], 2, stats::sd)
    cat("   names identical to the reference:",
        identical(names(w), colnames(keep[["integer"]])), "\n")
    cat("   max |row - reference| / sd(column):",
        signif(max(abs(got - w) / sp), 4), "\n")
    cat("   per coefficient |diff| / sd(column):",
        paste(signif(abs(got - w) / sp, 3), collapse = " "), "\n")
  }
  cat("\n")
}

dd <- mk(901)
arm(dd, which(dd$y == 2L), FALSE, "seed 901, ungrouped, one interior row")

dg2 <- mkg(902)
th2 <- c(a = 3L, b = 2L, c = 2L)
arm(dg2, which(dg2$y == 2L), TRUE,
    "seed 902, grouped, one interior row (level-local loss)", th2)

## the worker's seed 904: category 2 present on exactly one row, in level a,
## so the level disappears GLOBALLY when that row goes
mkg904 <- function(seed = 904, n = 96) {
  set.seed(seed)
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- stats::rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.0), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
  y <- vapply(seq_len(n), function(i) {
    1L + sum(stats::runif(1) >
               stats::plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
  }, 1L)
  mid <- which(y == 2L)
  y[mid[-1L]] <- 1L
  data.frame(x = x, g = g, y = y)
}
dg4 <- mkg904(904)
cat("seed 904 interior rows:", which(dg4$y == 2L),
    " in level(s)", as.character(dg4$g[dg4$y == 2L]), "\n")
arm(dg4, which(dg4$y == 2L)[1L], TRUE,
    "seed 904, grouped, one interior row globally", th2)
cat("DONE rev-20\n")
