## RE-CHECK item 4: the four influence() failure cases, and item 5, the
## control that a NON-threshold grouping factor deletion is untouched.
## Every case is run on both arms, and warnings are captured with their
## text rather than suppressed, because the count in the message is the
## thing under review.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, "\n\n")
options(width = 160)

run <- function(tag, expr) {
  cat("==", tag, "\n")
  ws <- character(0)
  out <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  cat("   warnings:", length(ws), "\n")
  for (w in ws) cat("     W:", gsub("\n", " ", w), "\n")
  if (inherits(out, "condition")) {
    cat("   ERROR (", class(out)[1L], "):\n     ",
        gsub("\n", " ", conditionMessage(out)), "\n\n")
    return(invisible(NULL))
  }
  if (inherits(out, "frmtmb_influence")) {
    m <- out$fixed
    filled <- which(apply(m, 1, function(r) !all(is.na(r))))
    cat("   table", paste(dim(m), collapse = "x"),
        " NA cells =", sum(is.na(m)),
        " rows entirely NA =", sum(apply(m, 1, function(r) all(is.na(r)))),
        " rows with any value =", length(filled), "\n")
    if (length(filled) && length(filled) <= 5L) {
      for (i in filled) {
        cat("     row", rownames(m)[i], ":",
            paste(signif(m[i, ], 8), collapse = "  "), "\n")
      }
    }
    cd <- suppressWarnings(cooks.distance(out))
    cat("   cooks NA =", sum(is.na(cd)), "of", length(cd), "\n\n")
  } else {
    print(out); cat("\n")
  }
  invisible(out)
}

## ---- case 1: data = whose response reaches a HIGHER category --------
set.seed(1101)
n <- 50
x <- stats::rnorm(n)
cp <- cbind(stats::plogis(-0.7 - 0.5 * x), stats::plogis(0.6 - 0.5 * x))
dd <- data.frame(x = x, y = 1L + rowSums(stats::runif(n) > cp))
fit <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = dd))
cat("fit n_tau =", length(fit$estimates[["tau_raw"]]),
    " table(y) =", paste(table(dd$y), collapse = "/"), "\n\n")
dhi <- dd; dhi$y[c(1L, 2L)] <- 4L
run("case 1: data = with TWO rows in a category the fit never saw",
    influence(fit, data = dhi, force = TRUE))

## ---- case 4 (partial): ONE row in an unseen category ---------------
dhi1 <- dd; dhi1$y[7L] <- 4L
run("case 4: data = with ONE row in an unseen category (row 7)",
    influence(fit, data = dhi1, force = TRUE))

## the ordered-factor spelling of the same thing, which goes through
## thres_pin_recode() rather than the count guard
ddf <- dd; ddf$y <- factor(ddf$y, levels = 1:3, ordered = TRUE)
fitf <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = ddf))
dhf <- ddf
dhf$y <- factor(c(as.character(ddf$y)[-7L][1:6], "4",
                  as.character(ddf$y)[-(1:7)]),
                levels = 1:4, ordered = TRUE)
dhf$y <- factor(as.character(ddf$y), levels = 1:4, ordered = TRUE)
dhf$y[7L] <- "4"
run("case 4b: the same as an ordered factor (recode route)",
    influence(fitf, data = dhf, force = TRUE))

## ---- case 2: data = with a thres(gr = ) level the fit never saw -----
set.seed(1103)
k <- 96
g <- factor(rep(c("a", "b", "c"), length.out = k))
xg <- stats::rnorm(k)
tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
yg <- vapply(seq_len(k), function(i) {
  1L + sum(stats::runif(1) >
             stats::plogis(tau[[as.character(g[i])]] - 0.5 * xg[i]))
}, 1L)
dg <- data.frame(x = xg, g = g, y = yg)
fg <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x), family = cumulative(),
                          data = dg))
cat("grouped fit nthres =",
    paste(fg$spec$responses$y$family[["thres"]][["nthres"]], collapse = "/"),
    "\n\n")
dnew <- dg
dnew$g <- as.character(dnew$g); dnew$g[1:8] <- "d"
dnew$g <- factor(dnew$g); dnew$y[1:8] <- pmin(dnew$y[1:8], 3L)
run("case 2: data = with a FOURTH thres(gr = ) level 'd'",
    influence(fg, data = dnew, force = TRUE))

## ---- case 3: groups = deleting a whole thres(gr = ) level ----------
dgr <- dg
fg2 <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x + (1 | g)),
                           family = cumulative(), data = dgr))
run("case 3: influence(groups = 'g') where g IS the threshold factor",
    influence(fg2, groups = "g"))

## ---- item 5: the control, a grouping factor that is NOT the
## threshold factor -------------------------------------------------
set.seed(1201)
m <- 120
gg <- factor(rep(c("a", "b", "c"), length.out = m))
id <- factor(rep(1:8, length.out = m))
xi <- stats::rnorm(m)
yi <- vapply(seq_len(m), function(i) {
  1L + sum(stats::runif(1) >
             stats::plogis(tau[[as.character(gg[i])]] - 0.5 * xi[i]))
}, 1L)
di <- data.frame(x = xi, g = gg, id = id, y = yi)
cat("control: per-level max =", paste(tapply(di$y, di$g, max),
                                      collapse = "/"), "\n")
cat("control: smallest rows any g level keeps over the 8 id deletions =",
    min(vapply(levels(di$id), function(l) {
      min(table(di$g[di$id != l]))
    }, 1L)), "\n\n")
fi <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x + (1 | id)),
                          family = cumulative(), data = di))
run("item 5 control: influence(groups = 'id'), id is NOT the thres factor",
    influence(fi, groups = "id"))
## and the same model's observation-wise deletion, nothing lost
run("item 5 control: observation-wise on the same model",
    influence(fi, force = TRUE))
## and the ungrouped control at seed 1101 with data untouched
run("item 5 control: seed 1101 with no doctored data",
    influence(fit, force = TRUE))
cat("DONE rev-22\n")
