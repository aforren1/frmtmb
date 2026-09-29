## RE-CHECK item 6: is the new error the right shape for a NON-ordinal fit
## whose data2 is lost, or should influence() refuse earlier? Measured:
## what a failed unit costs, and whether a one-call pre-flight on the FULL
## data would have caught each of the four failure cases.
## Also the inertness case rev-09 does not have: an ORDERED-FACTOR ordinal
## response with every category well populated, so no deletion empties one.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, "\n\n")
options(width = 150, digits = 17)

## ---- inertness with an ordered factor that loses nothing ------------
set.seed(2301)
n <- 80
x <- stats::rnorm(n)
cp <- cbind(stats::plogis(-0.9 - 0.5 * x), stats::plogis(0.1 - 0.5 * x),
            stats::plogis(1.1 - 0.5 * x))
y <- 1L + rowSums(stats::runif(n) > cp)
dd <- data.frame(x = x, y = y)
cat("table(y) =", paste(table(dd$y), collapse = "/"),
    " (every category has many rows, so no deletion empties one)\n")
tabs <- list()
for (how in c("integer", "ordered", "character")) {
  d <- dd
  d$y <- switch(how, integer = dd$y,
                ordered = factor(dd$y, levels = 1:4, ordered = TRUE),
                character = as.character(dd$y))
  fit <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = d))
  nw <- 0L
  inf <- withCallingHandlers(suppressWarnings(influence(fit, force = TRUE)),
                             warning = function(w) {
                               nw <<- nw + 1L; invokeRestart("muffleWarning")
                             })
  tabs[[how]] <- inf$fixed
  cat("  ", how, ": NA =", sum(is.na(inf$fixed)), " warnings =", nw,
      " sum =", format(sum(inf$fixed), digits = 17), "\n")
}
cat("   identical(ordered, integer):", identical(tabs$ordered, tabs$integer),
    "  identical(character, integer):",
    identical(tabs$character, tabs$integer), "\n\n")

## ---- would a pre-flight on the FULL data catch each failure case? ---
cat("#### a single assemble_frame() on the FULL data, per failure case\n")
probe <- function(tag, spec, data, pin, data2 = list()) {
  out <- tryCatch({
    if (is.null(pin)) {
      frmtmb:::assemble_frame(spec, data, data2 = data2)
    } else {
      frmtmb:::assemble_frame(spec, data, data2 = data2, thres_pin = pin)
    }
    "assembles"
  }, error = function(e) paste0("REFUSES: ",
                               substr(gsub("\n", " ",
                                           conditionMessage(e)), 1, 90)))
  cat("  ", tag, "->", out, "\n")
}
has_pin <- exists("thres_pin_of_fit", envir = asNamespace("frmtmb"),
                  inherits = FALSE)
set.seed(1101)
m <- 50
xm <- stats::rnorm(m)
cpm <- cbind(stats::plogis(-0.7 - 0.5 * xm), stats::plogis(0.6 - 0.5 * xm))
dm <- data.frame(x = xm, y = 1L + rowSums(stats::runif(m) > cpm))
f1 <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = dm))
p1 <- if (has_pin) frmtmb:::thres_pin_of_fit(f1) else NULL
dhi <- dm; dhi$y[c(1L, 2L)] <- 4L
probe("case 1 (data = , two unseen-category rows)", f1$spec, dhi, p1)
dh1 <- dm; dh1$y[7L] <- 4L
probe("case 4 (data = , one unseen-category row)", f1$spec, dh1, p1)

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
pg <- if (has_pin) frmtmb:::thres_pin_of_fit(fg) else NULL
dnew <- dg; dnew$g <- as.character(dnew$g); dnew$g[1:8] <- "d"
dnew$g <- factor(dnew$g); dnew$y[1:8] <- pmin(dnew$y[1:8], 3L)
probe("case 2 (data = , a fourth thres(gr = ) level)", fg$spec, dnew, pg)
probe("case 3 (groups = 'g'), the FULL data", fg$spec, dg, pg)
cat("   case 3 is the one a pre-flight on the full data cannot see: the\n",
    "  full data assembles and only the per-level subsets do not.\n\n")

## ---- what does a failed unit cost? ---------------------------------
cat("#### the cost of the failures the counter now reports\n")
t0 <- proc.time()
r <- tryCatch(influence(f1, data = dhi, force = TRUE),
              error = function(e) e)
t1 <- proc.time()
cat("   50 failing units (they die at frame assembly, before any fit):",
    sprintf("%.3f", (t1 - t0)[["elapsed"]]), "s\n")
t0 <- proc.time()
ok <- suppressWarnings(influence(f1, force = TRUE))
t1 <- proc.time()
cat("   50 SUCCEEDING units, the same model, for scale:",
    sprintf("%.3f", (t1 - t0)[["elapsed"]]), "s\n")
cat("   so the whole set of failures costs a fraction of one ordinary\n",
    "  influence() call; a pre-flight would save that, not more.\n")
cat("DONE rev-23\n")
