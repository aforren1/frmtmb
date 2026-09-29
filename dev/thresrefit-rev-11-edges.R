## REVIEW, edges the findings assert but do not show, and edges the pin
## could break: influence(data = ) with a HIGHER top category, a
## multivariate model with two ordinal responses, influence(groups = ) on
## a thres(gr = ) model, and a pin met by a grouping level the fit never
## saw (the refusal thres_pin_apply() writes).
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, " ver =",
    as.character(utils::packageVersion("frmtmb")), "\n\n")
options(width = 150)

say <- function(tag, expr) {
  cat("==", tag, "\n")
  out <- tryCatch(expr, error = function(e) e, warning = function(w) w)
  if (inherits(out, "condition")) {
    cat("   ", class(out)[1L], ":",
        gsub("\n", " ", conditionMessage(out)), "\n\n")
  } else {
    cat("    returned without error\n"); print(out); cat("\n")
  }
  invisible(out)
}

## ---- 1. influence(data = ) whose response reaches a HIGHER category --
set.seed(1101)
n <- 50
x <- stats::rnorm(n)
cp <- cbind(stats::plogis(-0.7 - 0.5 * x), stats::plogis(0.6 - 0.5 * x))
dd <- data.frame(x = x, y = 1L + rowSums(stats::runif(n) > cp))
cat("fitted table(y) =", paste(table(dd$y), collapse = "/"), "\n")
fit <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = dd))
cat("fitted n_tau =", length(fit$estimates[["tau_raw"]]), "\n")
dhi <- dd
dhi$y[c(1L, 2L)] <- 4L      # two rows in a category the fit never saw
say("influence(data = ) whose response reaches category 4",
    {
      inf <- suppressWarnings(influence(fit, data = dhi, force = TRUE))
      c(ncol = ncol(inf$fixed), na = sum(is.na(inf$fixed)),
        rows_all_na = sum(apply(inf$fixed, 1, function(r) all(is.na(r)))))
    })

## ---- 2. a multivariate model with two ordinal responses -------------
set.seed(1102)
m <- 60
xm <- stats::rnorm(m)
cpa <- cbind(stats::plogis(-0.6 - 0.5 * xm), stats::plogis(0.5 - 0.5 * xm),
             stats::plogis(2.4 - 0.5 * xm))
ya <- 1L + rowSums(stats::runif(m) > cpa)
ta <- which(ya == 4L); if (length(ta) > 1L) ya[ta[-1L]] <- 3L
cpb <- cbind(stats::plogis(-0.4 - 0.3 * xm), stats::plogis(0.8 - 0.3 * xm))
yb <- 1L + rowSums(stats::runif(m) > cpb)
dm <- data.frame(x = xm, ya = ya, yb = yb)
cat("mv tables:", paste(table(dm$ya), collapse = "/"), " and ",
    paste(table(dm$yb), collapse = "/"), "\n")
say("two ordinal responses, influence()",
    {
      fm <- suppressWarnings(frm(bf(ya ~ x) + bf(yb ~ x),
                                family = cumulative(), data = dm))
      if (exists("thres_pin_of_fit", envir = asNamespace("frmtmb"),
                 inherits = FALSE)) {
        p <- frmtmb:::thres_pin_of_fit(fm)
        cat("    pin names:", paste(names(p), collapse = ","),
            " nthres:", paste(vapply(p, function(z) z$nthres, 1L),
                              collapse = "/"), "\n")
      }
      inf <- suppressWarnings(influence(fm, force = TRUE))
      c(ncol = ncol(inf$fixed), na = sum(is.na(inf$fixed)),
        na_top = sum(is.na(inf$fixed[which(dm$ya == 4L), ])))
    })

## ---- 3. influence(groups = ) on a thres(gr = ) model -----------------
set.seed(1103)
k <- 96
g <- factor(rep(c("a", "b", "c"), length.out = k))
xg <- stats::rnorm(k)
tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
yg <- vapply(seq_len(k), function(i) {
  1L + sum(stats::runif(1) >
             stats::plogis(tau[[as.character(g[i])]] - 0.5 * xg[i]))
}, 1L)
dg <- data.frame(x = xg, g = g, y = yg, sub = factor(rep(1:12, 8)))
say("influence(groups = 'sub') on a thres(gr = g) model",
    {
      fg <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x + (1 | sub)),
                                family = cumulative(), data = dg))
      inf <- suppressWarnings(influence(fg, groups = "sub"))
      c(nrow = nrow(inf$fixed), ncol = ncol(inf$fixed),
        na = sum(is.na(inf$fixed)))
    })
say("influence(groups = 'g') on a thres(gr = g) model: a level goes",
    {
      fg2 <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x + (1 | g)),
                                 family = cumulative(), data = dg))
      inf <- suppressWarnings(influence(fg2, groups = "g"))
      c(nrow = nrow(inf$fixed), ncol = ncol(inf$fixed),
        na = sum(is.na(inf$fixed)))
    })

## ---- 4. the refusal thres_pin_apply() writes, reached by construction
if (exists("thres_pin_apply", envir = asNamespace("frmtmb"),
           inherits = FALSE)) {
  gi <- structure(c(1, 1, 2, 2), thres_levels = c("a", "zz"))
  say("thres_pin_apply() met by a grouping level the fit never saw",
      frmtmb:::thres_pin_apply(
        list(y = list(grouped = TRUE, nthres = c(3L, 2L),
                      groups = c("a", "b"), levels_y = NULL)),
        list(resp_name = "y"), list(thres_gr = gi)))
  ## and the same through influence(data = ) with a new level
  say("influence(data = ) with a thres(gr = ) level the fit never saw",
      {
        fg <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                                  family = cumulative(), data = dg))
        dnew <- dg
        dnew$g <- as.character(dnew$g)
        dnew$g[1:8] <- "d"
        dnew$g <- factor(dnew$g)
        inf <- suppressWarnings(influence(fg, data = dnew, force = TRUE))
        c(ncol = ncol(inf$fixed), na = sum(is.na(inf$fixed)))
      })
}
cat("DONE rev-11\n")
