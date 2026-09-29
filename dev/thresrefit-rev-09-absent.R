## REVIEW claim 6, the guard's ABSENT case: a family with nothing to pin
## must reach thres_pin_of_fit() and get NULL, and influence() must be
## bit-identical to the base build. Printed at 17 digits and as a
## checksum, because "%.6f" renders 1.8e-11 as 0.000000.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, " ver =",
    as.character(utils::packageVersion("frmtmb")), "\n\n")
options(width = 140, digits = 17)

has_pin <- exists("thres_pin_of_fit", envir = asNamespace("frmtmb"),
                  inherits = FALSE)
cat("thres_pin_of_fit present:", has_pin, "\n\n")

report <- function(tag, fit) {
  if (has_pin) {
    p <- frmtmb:::thres_pin_of_fit(fit)
    cat("  thres_pin_of_fit ->",
        if (is.null(p)) "NULL" else
          paste0("list of ", length(p), ": ",
                 paste(vapply(p, function(z) paste0("grouped=", z$grouped,
                                                    " nthres=",
                                                    paste(z$nthres,
                                                          collapse = "/")),
                              ""), collapse = "; ")), "\n")
  }
  inf <- try(suppressWarnings(influence(fit, force = TRUE)), silent = TRUE)
  if (inherits(inf, "try-error")) {
    cat("  influence() REFUSED:",
        gsub("\n", " ", conditionMessage(attr(inf, "condition"))), "\n")
    return(invisible(NULL))
  }
  cat("  influence dim:", paste(dim(inf$fixed), collapse = "x"),
      " NA cells:", sum(is.na(inf$fixed)), "\n")
  cat("  colnames:", paste(colnames(inf$fixed), collapse = ","), "\n")
  ## a digest of the whole table at full precision: one changed bit shows
  cat("  md5(fixed):",
      digest_or_sum(inf$fixed), "\n")
  cat("  sum(fixed):", format(sum(inf$fixed), digits = 17), "\n")
  cat("  cooks[1:4]:",
      paste(format(utils::head(suppressWarnings(cooks.distance(inf)), 4),
                   digits = 17), collapse = " "), "\n")
  invisible(NULL)
}
digest_or_sum <- function(m) {
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(round(m, 12))
  } else {
    format(sum(abs(m)), digits = 17)
  }
}

set.seed(601)
n <- 40
d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:8, 5)))
d$yg <- stats::rnorm(n, 1 + d$x)
d$yb <- stats::rbinom(n, 1, stats::plogis(d$x))
d$yp <- stats::rpois(n, exp(0.4 + 0.3 * d$x))
d$yc <- factor(sample(c("a", "b", "c"), n, TRUE))

cat("== gaussian, no random effect\n")
report("gauss", frm(bf(yg ~ x), family = gaussian(), data = d))
cat("== gaussian with (1 | g)\n")
report("gaussre", suppressWarnings(frm(bf(yg ~ x + (1 | g)),
                                       family = gaussian(), data = d)))
cat("== binomial\n")
report("binom", suppressWarnings(frm(bf(yb ~ x), family = bernoulli(),
                                     data = d)))
cat("== poisson\n")
report("pois", suppressWarnings(frm(bf(yp ~ x), family = poisson(),
                                    data = d)))
cat("== categorical\n")
report("cat", suppressWarnings(frm(bf(yc ~ x), family = categorical(),
                                   data = d)))
cat("== ordinal WITHOUT a lost category (the pin present and inert)\n")
set.seed(602)
do <- data.frame(x = stats::rnorm(60))
cpo <- cbind(stats::plogis(-0.7 - 0.5 * do$x),
             stats::plogis(0.6 - 0.5 * do$x),
             stats::plogis(1.4 - 0.5 * do$x))
do$y <- 1L + rowSums(stats::runif(60) > cpo)
cat("   table(y) =", paste(table(do$y), collapse = "/"), "\n")
report("ord_full", suppressWarnings(frm(bf(y ~ x), family = cumulative(),
                                        data = do)))
cat("== ordinal with thres(3) written by hand, no lost category\n")
report("ord_pin", suppressWarnings(frm(bf(y | thres(3) ~ x),
                                       family = cumulative(), data = do)))
cat("== ordinal grouped, no lost category\n")
set.seed(603)
dgg <- data.frame(g = factor(rep(c("a", "b"), 45)), x = stats::rnorm(90))
tt <- list(a = c(-0.6, 0.6), b = c(-0.4, 0.8))
dgg$y <- vapply(seq_len(90), function(i) {
  1L + sum(stats::runif(1) >
             stats::plogis(tt[[as.character(dgg$g[i])]] - 0.5 * dgg$x[i]))
}, 1L)
cat("   per-level max =", paste(tapply(dgg$y, dgg$g, max), collapse = "/"),
    "\n")
report("ord_gr", suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                                      family = cumulative(), data = dgg)))
cat("DONE rev-09\n")
