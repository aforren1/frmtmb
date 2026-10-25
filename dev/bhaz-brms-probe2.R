# brms 2.23.0's bhaz() with knots = and Boundary.knots =. brms checks
# every variable an addition term names against data and data2, then
# evaluates the flags in the formula's environment, so literal vectors
# are the spelling that works. Same data as dev/bhaz-brms-probe.R.
#   Rscript dev/bhaz-brms-probe2.R > dev/bhaz-log/brms-probe2.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(11); n <- 200; x <- rnorm(n)
g <- factor(sample(c("a", "b", "c"), n, TRUE))
tt <- rexp(n, exp(-0.5 + 0.7 * x)); ct <- rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                 g = g)
ns <- asNamespace("frmtmb.survival")

basis <- function(lab, fo, data2 = NULL) {
  cat("\n=====", lab, "=====\n")
  r <- tryCatch(withCallingHandlers(
    standata(fo, data = dd, family = cox(), data2 = data2),
    warning = function(w) {
      cat("warning:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(r, "error")) {
    cat("error:", conditionMessage(r), "\n")
    return(invisible(NULL))
  }
  at <- attributes(r$Zbhaz)
  cat("Kbhaz", r$Kbhaz, " degree", at$degree, " intercept", at$intercept,
      "\n knots", format(at$knots, digits = 17),
      "\n Boundary.knots", format(at$Boundary.knots, digits = 17), "\n")
  cat(" rows of Zbhaz that are all zero:", sum(rowSums(r$Zbhaz) == 0),
      "  rows with a negative entry:", sum(apply(r$Zbhaz, 1, min) < 0),
      "\n rows of Zcbhaz all zero:", sum(rowSums(r$Zcbhaz) == 0),
      "  max row sum of Zcbhaz", format(max(rowSums(r$Zcbhaz)),
                                         digits = 17), "\n")
  invisible(r)
}

k1 <- basis("bhaz(knots = c(0.5, 1.1), Boundary.knots = c(0, 16))",
            bf(time | cens(cens) + bhaz(knots = c(0.5, 1.1),
                                        Boundary.knots = c(0, 16)) ~ x))
k2 <- basis("bhaz(knots = c(0.5, 1.1))",
            bf(time | cens(cens) + bhaz(knots = c(0.5, 1.1)) ~ x))
k3 <- basis("bhaz(df = 8, knots = c(0.5, 1.1))",
            bf(time | cens(cens) + bhaz(df = 8, knots = c(0.5, 1.1)) ~ x))
k4 <- basis("bhaz(Boundary.knots = c(0, 16))",
            bf(time | cens(cens) + bhaz(Boundary.knots = c(0, 16)) ~ x))
k5 <- basis("bhaz(df = 6, Boundary.knots = c(0, 16))",
            bf(time | cens(cens) + bhaz(df = 6, Boundary.knots = c(0, 16))
               ~ x))
k6 <- basis("bhaz(Boundary.knots = c(0.5, 10)): inside the data",
            bf(time | cens(cens) + bhaz(Boundary.knots = c(0.5, 10)) ~ x))
k7 <- basis("bhaz(knots = c(0.5, 20)): a knot outside the data",
            bf(time | cens(cens) + bhaz(knots = c(0.5, 20)) ~ x))
k8 <- basis("bhaz(knots = numeric(0), df = 4)",
            bf(time | cens(cens) + bhaz(knots = numeric(0), df = 4) ~ x))
kn <- c(0.5, 1.1)
k9 <- basis("bhaz(knots = kn) with kn in data2 and the env",
            bf(time | cens(cens) + bhaz(knots = kn) ~ x),
            data2 = list(kn = kn))
k10 <- basis("bhaz(df = 4.5)",
             bf(time | cens(cens) + bhaz(df = 4.5) ~ x))
k11 <- basis("bhaz(df = 3)", bf(time | cens(cens) + bhaz(df = 3) ~ x))
k12 <- basis("bhaz(df = 2)", bf(time | cens(cens) + bhaz(df = 2) ~ x))
k13 <- basis("bhaz(foo = 1): an argument splines2 does not take",
             bf(time | cens(cens) + bhaz(foo = 1) ~ x))
k14 <- basis("bhaz(gr = g, knots = c(0.5, 1.1))",
             bf(time | cens(cens) + bhaz(gr = g, knots = c(0.5, 1.1)) ~ x))
k15 <- basis("bhaz(df = 6) + bhaz(df = 7): twice",
             bf(time | cens(cens) + bhaz(df = 6) + bhaz(df = 7) ~ x))

# frmtmb.survival 0.1.0's basis, on a spec built from brms's own knots
cmp <- function(lab, r, degree = 3L, intercept = TRUE) {
  if (is.null(r)) return(invisible(NULL))
  at <- attributes(r$Zbhaz)
  sp <- list(internal = as.numeric(at$knots),
             boundary = as.numeric(at$Boundary.knots), degree = degree,
             intercept = intercept, df = r$Kbhaz)
  y <- dd$time
  inside <- y >= sp$boundary[1] & y <= sp$boundary[2]
  M <- ns$mspline_design(y, sp)
  I <- ns$ispline_design(y, sp)
  bz <- matrix(unclass(r$Zbhaz), nrow(r$Zbhaz))
  cz <- matrix(unclass(r$Zcbhaz), nrow(r$Zcbhaz))
  cat(sprintf("%-55s rows inside %d of %d  max|dM|/max|M| %.3e  max|dI|/max|I| %.3e\n",
              lab, sum(inside), length(y),
              max(abs(M[inside, ] - bz[inside, ])) / max(abs(bz)),
              max(abs(I[inside, ] - cz[inside, ])) / max(abs(cz))))
}
cat("\n===== frmtmb.survival 0.1.0's basis on brms's knots =====\n")
cmp("knots = c(0.5, 1.1), Boundary.knots = c(0, 16)", k1)
cmp("knots = c(0.5, 1.1)", k2)
cmp("df = 8, knots = c(0.5, 1.1)", k3)
cmp("Boundary.knots = c(0, 16)", k4)
cmp("df = 6, Boundary.knots = c(0, 16)", k5)
cmp("Boundary.knots = c(0.5, 10)", k6)
