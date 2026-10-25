# Smoke test of the bhaz() spellings on the lane build.
.libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb.survival))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
set.seed(11); n <- 200; x <- rnorm(n)
g <- factor(sample(c("a", "b", "c"), n, TRUE))
tt <- rexp(n, exp(-0.5 + 0.7 * x)); ct <- rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                 g = g)
try1 <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-50s ERROR: %s\n", lab, conditionMessage(r)))
  } else if (inherits(r, "frmtmb_fit")) {
    cat(sprintf("%-50s logLik %.12f  b_x %.10f\n", lab,
                as.numeric(logLik(r)), fixef(r)["x", 1]))
  } else {
    cat(sprintf("%-50s OK (%s)\n", lab, class(r)[1]))
  }
  invisible(r)
}
a <- try1("cox(df = 6)", frm(bf(time | cens(cens) ~ x),
                             family = cox(df = 6), data = dd))
b <- try1("bhaz(df = 6)", frm(bf(time | cens(cens) + bhaz(df = 6) ~ x),
                              family = cox(), data = dd))
cat("identical logLik:", identical(logLik(a), logLik(b)),
    " identical estimates:", identical(a$estimates, b$estimates), "\n")
try1("bhaz(df = 6) + cox(df = 6)",
     frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = cox(df = 6),
         data = dd))
try1("bhaz(df = 6) + cox(df = 5)",
     frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = cox(df = 5),
         data = dd))
try1("bhaz(knots, Boundary.knots)",
     frm(bf(time | cens(cens) + bhaz(knots = c(0.5, 1.1),
                                     Boundary.knots = c(0, 16)) ~ x),
         family = cox(), data = dd))
try1("cox(knots, Boundary.knots)",
     frm(bf(time | cens(cens) ~ x),
         family = cox(knots = c(0.5, 1.1), Boundary.knots = c(0, 16)),
         data = dd))
kn <- c(0.5, 1.1)
try1("bhaz(knots = kn) from the env",
     frm(bf(time | cens(cens) + bhaz(knots = kn) ~ x), family = cox(),
         data = dd))
try1("bhaz(df = 8, knots)",
     frm(bf(time | cens(cens) + bhaz(df = 8, knots = c(0.5, 1.1)) ~ x),
         family = cox(), data = dd))
try1("bhaz(knots) + cox(df = 4)",
     frm(bf(time | cens(cens) + bhaz(knots = c(0.5, 1.1)) ~ x),
         family = cox(df = 4), data = dd))
try1("bhaz(gr = g)", frm(bf(time | cens(cens) + bhaz(gr = g) ~ x),
                         family = cox(), data = dd))
try1("bhaz(g)", frm(bf(time | cens(cens) + bhaz(g) ~ x),
                    family = cox(), data = dd))
try1("bhaz(gr = NA, df = 6)",
     frm(bf(time | cens(cens) + bhaz(gr = NA, df = 6) ~ x),
         family = cox(), data = dd))
try1("bhaz(foo = 1)", frm(bf(time | cens(cens) + bhaz(foo = 1) ~ x),
                          family = cox(), data = dd))
try1("bhaz(Boundary.knots = c(0.5, 10))",
     frm(bf(time | cens(cens) + bhaz(Boundary.knots = c(0.5, 10)) ~ x),
         family = cox(), data = dd))
try1("bhaz(knots = c(0.5, 20))",
     frm(bf(time | cens(cens) + bhaz(knots = c(0.5, 20)) ~ x),
         family = cox(), data = dd))
try1("gaussian + bhaz(df = 6)",
     frm(bf(time | bhaz(df = 6) ~ x), family = gaussian(), data = dd))
try1("bhaz twice",
     frm(bf(time | cens(cens) + bhaz(df = 6) + bhaz(df = 6) ~ x),
         family = cox(), data = dd))
try1("cox(bhaz = list(df = 6))", cox(bhaz = list(df = 6)))
try1("bhaz(df = 3)", frm(bf(time | cens(cens) + bhaz(df = 3) ~ x),
                         family = cox(), data = dd))
try1("bhaz(df = nope)", frm(bf(time | cens(cens) + bhaz(df = nope) ~ x),
                            family = cox(), data = dd))
try1("family = 'cox' + bhaz(df = 6)",
     frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = "cox",
         data = dd))
try1("brmsfamily('cox') + bhaz(df = 6)",
     frm(bf(time | cens(cens) + bhaz(df = 6) ~ x),
         family = brmsfamily("cox"), data = dd))
# update(), predict() on new data and a refit reparse the formula
u <- try1("update(b, . ~ . + g)", update(b, . ~ . + g))
p <- try1("frm_linpred(b, newdata)",
          frm_linpred(b, newdata = dd[1:5, ], type = "link"))
print(head(p))
cat("Kbhaz of b:", length(cox_baseline(b)), "\n")
print(family(b))
