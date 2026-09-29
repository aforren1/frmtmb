## RE-CHECK item 8: are the two new "cannot be refit" statements true of
## the CODE, and are they no wider than the code?
##  (a) a fit with thres(K) pinned by hand declares categories it never
##      saw. Does influence(data = ) holding one of them refuse, as the
##      sentence says, or succeed, as the model allows?
##  (b) the Details say "top or interior". What about the BOTTOM category?
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, "\n\n")
options(width = 150, digits = 17)

run <- function(tag, expr) {
  cat("==", tag, "\n")
  ws <- character(0)
  out <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                             warning = function(w) {
                               ws <<- c(ws, conditionMessage(w))
                               invokeRestart("muffleWarning")
                             })
  for (w in ws) cat("   W:", gsub("\n", " ", w), "\n")
  if (inherits(out, "condition")) {
    cat("   ERROR:", gsub("\n", " ", conditionMessage(out)), "\n\n")
    return(invisible(NULL))
  }
  if (inherits(out, "frmtmb_influence")) {
    cat("   table", paste(dim(out$fixed), collapse = "x"),
        " NA cells =", sum(is.na(out$fixed)),
        " rows entirely NA =",
        sum(apply(out$fixed, 1, function(r) all(is.na(r)))),
        " warnings =", length(ws), "\n\n")
  } else {
    print(out); cat("\n")
  }
  invisible(out)
}

## ---- (a) a hand-pinned count that exceeds the observed categories ----
set.seed(2401)
n <- 60
x <- stats::rnorm(n)
cp <- cbind(stats::plogis(-0.7 - 0.5 * x), stats::plogis(0.8 - 0.5 * x))
dd <- data.frame(x = x, y = 1L + rowSums(stats::runif(n) > cp))
cat("data table(y) =", paste(table(dd$y), collapse = "/"),
    " (max 3), and the model pins thres(4), so 5 categories\n")
f4 <- try(suppressWarnings(frm(bf(y | thres(4) ~ x), family = cumulative(),
                               data = dd)), silent = TRUE)
if (inherits(f4, "try-error")) {
  cat("   frm(thres(4)) REFUSED:",
      gsub("\n", " ", conditionMessage(attr(f4, "condition"))), "\n")
} else {
  cat("   fitted n_tau =", length(f4$estimates[["tau_raw"]]), "\n")
  d4 <- dd; d4$y[1L] <- 4L   # a category the FIT never saw but the MODEL has
  run(paste0("(a) influence(data = ) with one row in category 4, a category ",
             "the fit never saw but thres(4) declares"),
      influence(f4, data = d4, force = TRUE))
  d5 <- dd; d5$y[1L] <- 6L   # beyond the model as well
  run("(a2) the same with category 6, beyond thres(4) too",
      influence(f4, data = d5, force = TRUE))
  ## the ordered-factor spelling cannot even be fitted, for the record
  do <- dd; do$y <- factor(do$y, levels = 1:5, ordered = TRUE)
  fo <- try(suppressWarnings(frm(bf(y | thres(4) ~ x),
                                 family = cumulative(), data = do)),
            silent = TRUE)
  cat("== (a3) the same model with an ordered factor of 5 declared levels:",
      if (inherits(fo, "try-error"))
        paste("REFUSED:",
              substr(gsub("\n", " ",
                          conditionMessage(attr(fo, "condition"))), 1, 110))
      else paste("fitted, n_tau =", length(fo$estimates[["tau_raw"]])), "\n\n")
}

## ---- (b) the BOTTOM category, which the Details do not name ----------
set.seed(2402)
m <- 60
xm <- stats::rnorm(m)
cpm <- cbind(stats::plogis(-2.6 - 0.5 * xm), stats::plogis(-0.4 - 0.5 * xm),
             stats::plogis(0.9 - 0.5 * xm))
ym <- 1L + rowSums(stats::runif(m) > cpm)
lo <- which(ym == 1L)
if (length(lo) < 2L) stop("seed gives too few category-1 rows")
ym[lo[-1L]] <- 2L        # exactly one row left in the BOTTOM category
db <- data.frame(x = xm, y = ym)
ilo <- which(db$y == 1L)
cat("bottom-category construction: table(y) =",
    paste(table(db$y), collapse = "/"), " single bottom row =", ilo, "\n")
tb <- list()
for (how in c("integer", "ordered", "character")) {
  d <- db
  d$y <- switch(how, integer = db$y,
                ordered = factor(db$y, levels = 1:4, ordered = TRUE),
                character = as.character(db$y))
  fit <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = d))
  nw <- 0L
  inf <- withCallingHandlers(suppressWarnings(influence(fit, force = TRUE)),
                             warning = function(w) {
                               nw <<- nw + 1L; invokeRestart("muffleWarning")
                             })
  tb[[how]] <- inf$fixed
  cat("  ", how, ": n_tau =", length(fit$estimates[["tau_raw"]]),
      " NA =", sum(is.na(inf$fixed)), " warnings =", nw,
      " NA in the bottom row =", sum(is.na(inf$fixed[ilo, ])), "\n")
  cat("       bottom row:", paste(signif(inf$fixed[ilo, ], 9),
                                  collapse = "  "), "\n")
}
cat("   identical(ordered, integer):", identical(tb$ordered, tb$integer),
    "  identical(character, integer):",
    identical(tb$character, tb$integer), "\n")
ref <- suppressWarnings(frm(bf(y | thres(3) ~ x), family = cumulative(),
                           data = db[-ilo, , drop = FALSE]))
w <- frmtmb:::get_coef.frmtmb_fit(ref)
cat("   thres(3) reference:", paste(signif(w, 9), collapse = "  "), "\n")
cat("   names match:", identical(names(w), colnames(tb$integer)), "\n")
sp <- apply(tb$integer, 2, stats::sd)
cat("   per coefficient |row - reference| / sd(column):",
    paste(signif(abs(tb$integer[ilo, ] - w) / sp, 3), collapse = " "), "\n")
cat("DONE rev-24\n")
