# What brms 2.23.0 does with partial and unknown argument names in
# bhaz(). Same data as dev/bhaz-brms-probe.R (seed 11, n = 200, g drawn
# before the times; times 0.0044 to 14.42). No compile.
#   Rscript dev/bhaz-brms-partial.R > dev/bhaz-log/brms-partial.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(11); n <- 200; x <- rnorm(n)
g <- factor(sample(c("a", "b", "c"), n, TRUE))
tt <- rexp(n, exp(-0.5 + 0.7 * x)); ct <- rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                 g = g)
ref <- standata(bf(time | cens(cens) ~ x), data = dd, family = cox())
show <- function(term) {
  fo <- as.formula(paste0("time | cens(cens) + ", term, " ~ x"))
  r <- tryCatch(withCallingHandlers(
    standata(bf(fo), data = dd, family = cox()),
    warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-42s ERROR %s\n", term, conditionMessage(r)))
    return(invisible())
  }
  at <- attributes(r$Zbhaz)
  same <- identical(dim(r$Zbhaz), dim(ref$Zbhaz)) &&
    max(abs(unclass(r$Zbhaz) - unclass(ref$Zbhaz))) == 0
  cat(sprintf(paste0("%-42s K %d degree %d intercept %s knots [%s] ",
                     "boundary [%s] ngr %s same-as-default %s\n"),
              term, r$Kbhaz, at$degree, at$intercept,
              paste(format(at$knots, digits = 4), collapse = " "),
              paste(format(at$Boundary.knots, digits = 4),
                    collapse = " "),
              if (is.null(r$ngrbhaz)) "-" else r$ngrbhaz, same))
}
for (term in c("bhaz()", "bhaz(df = 6)", "bhaz(d = 6)", "bhaz(g = g)",
               "bhaz(intercept = FALSE)", "bhaz(int = FALSE)",
               "bhaz(inter = FALSE)", "bhaz(i = FALSE)",
               "bhaz(int = TRUE)", "bhaz(integral = TRUE)",
               "bhaz(Boundary.knots = c(0, 30))",
               "bhaz(Boundary = c(0, 30))", "bhaz(Bound = c(0, 30))",
               "bhaz(B = c(0, 30))", "bhaz(knots = c(0.5, 1.1))",
               "bhaz(kn = c(0.5, 1.1))", "bhaz(k = c(0.5, 1.1))",
               "bhaz(degree = 2)", "bhaz(deg = 2)", "bhaz(de = 2)",
               "bhaz(warn.outside = FALSE)", "bhaz(warn = FALSE)",
               "bhaz(w = FALSE)", "bhaz(per = TRUE)",
               "bhaz(periodic = TRUE)", "bhaz(derivs = 1)",
               "bhaz(foo = 1)", "bhaz(Knots = c(0.5, 1.1))",
               "bhaz(NA, 6)", "bhaz(NA, 6, c(0.5, 1.1))", "bhaz(6)",
               "bhaz(df = 6, d = 7)")) {
  show(term)
}
# Fix round 2: a full name and its shortened form together
for (term in c("bhaz(degree = 2, deg = 3)",
               "bhaz(knots = c(0.5, 1.1), kn = c(1, 2))",
               "bhaz(df = 6, df = 7)")) {
  show(term)
}
