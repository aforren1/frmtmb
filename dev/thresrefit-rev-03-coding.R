## REVIEW claim 3: with an ordered-factor response the model frame drops
## the unused TOP level. Check that influence() gives the SAME rows for
## integer, ordered factor, unordered factor and character codings, with
## and without thres(gr = g). Worker's seeds 501 and 502.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, " ver =",
    as.character(utils::packageVersion("frmtmb")), "\n\n")
options(width = 150)

refit_one_top <- function(seed = 501, n = 50) {
  set.seed(seed)
  x <- stats::rnorm(n)
  cp <- cbind(stats::plogis(-0.7 - 0.5 * x), stats::plogis(0.6 - 0.5 * x),
              stats::plogis(2.3 - 0.5 * x))
  y <- 1L + rowSums(stats::runif(n) > cp)
  top <- which(y == 4L)
  y[top[-1L]] <- 3L
  data.frame(x = x, y = y)
}
refit_grouped_top <- function(seed = 502, n = 96) {
  set.seed(seed)
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- stats::rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
  y <- vapply(seq_len(n), function(i) {
    1L + sum(stats::runif(1) >
               stats::plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
  }, 1L)
  ia <- which(g == "a" & y == 4L)
  y[ia[-1L]] <- 3L
  data.frame(x = x, g = g, y = y)
}

recode <- function(y, how) switch(how,
  integer = y,
  ordered = factor(y, levels = sort(unique(y)), ordered = TRUE),
  unordered = factor(y, levels = sort(unique(y))),
  character = as.character(y))

run <- function(dd, form, tag, irow) {
  cat("####", tag, "  deleted-unit row =", irow, "\n")
  keep <- list()
  for (how in c("integer", "ordered", "unordered", "character")) {
    d <- dd
    d$y <- recode(dd$y, how)
    fit <- try(suppressWarnings(frm(form, family = cumulative(), data = d)),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      cat("  ", how, ": frm() REFUSED:",
          substr(sub("\n.*", "", conditionMessage(attr(fit, "condition"))),
                 1, 90), "\n")
      next
    }
    inf <- try(suppressWarnings(influence(fit, force = TRUE)), silent = TRUE)
    if (inherits(inf, "try-error")) {
      cat("  ", how, ": influence() REFUSED\n"); next
    }
    keep[[how]] <- inf$fixed
    cat("  ", how, ": n_tau =", length(fit$estimates[["tau_raw"]]),
        " ncol =", ncol(inf$fixed), " NA cells =", sum(is.na(inf$fixed)),
        " NA in deleted row =", sum(is.na(inf$fixed[irow, ])), "\n")
    cat("       row:", paste(signif(inf$fixed[irow, ], 9), collapse = "  "),
        "\n")
    cat("       cooks[row] =", suppressWarnings(cooks.distance(inf))[irow],
        "\n")
  }
  if (length(keep) >= 2L) {
    base <- keep[["integer"]]
    sp <- apply(base, 2, stats::sd)
    for (how in setdiff(names(keep), "integer")) {
      m <- keep[[how]]
      if (!identical(dim(m), dim(base))) {
        cat("   dim mismatch integer vs", how, "\n"); next
      }
      d1 <- max(abs(m[irow, ] - base[irow, ]) / sp, na.rm = TRUE)
      dall <- max(abs(m - base) / rep(sp, each = nrow(base)), na.rm = TRUE)
      cat("   max |", how, "- integer| / sd(column): deleted row =",
          signif(d1, 4), "  whole table =", signif(dall, 4),
          "  identical() =", identical(m, base), "\n")
    }
  }
  cat("\n")
}

dd <- refit_one_top()
cat("ungrouped table(y) =", paste(table(dd$y), collapse = "/"), "\n")
run(dd, bf(y ~ x), "ungrouped, no thres()", which(dd$y == 4L))
run(dd, bf(y | thres(3) ~ x), "ungrouped, thres(3) written by hand",
    which(dd$y == 4L))

dg <- refit_grouped_top()
cat("grouped table(y) =", paste(table(dg$y), collapse = "/"),
    " per-level max =", paste(tapply(dg$y, dg$g, max), collapse = "/"), "\n")
run(dg, bf(y | thres(gr = g) ~ x), "grouped thres(gr = g)",
    which(dg$g == "a" & dg$y == 4L))
cat("DONE rev-03\n")
