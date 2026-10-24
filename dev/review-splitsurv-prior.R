# When is set_prior("dirichlet(2)", class = "sigma") refused: frmtmb
# base, frmtmb lane, brms 2.23.0. Usage: Rscript ... <base|lane>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-splitsurv-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
at <- function(lab, expr) {
  r <- tryCatch({ force(expr); "no error" },
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-40s %s\n", lab, substr(r, 1, 200)))
}
set.seed(1)
dd <- data.frame(y = rnorm(50), x = rnorm(50))
at("frmtmb set_prior(class=sigma)",
   p1 <- set_prior("dirichlet(2)", class = "sigma"))
at("frmtmb set_prior(class=Intercept,dpar=sigma)",
   p2 <- set_prior("dirichlet(2)", class = "Intercept", dpar = "sigma"))
at("frmtmb set_prior(class=sd)", set_prior("dirichlet(2)", class = "sd"))
at("frmtmb set_prior(class=foo)",
   p3 <- set_prior("dirichlet(2)", class = "foo"))
if (exists("p1")) {
  at("frmtmb validate_prior(p1)", validate_prior(p1, y ~ x, data = dd))
  at("frmtmb get_prior-free frm(p1)", frm(y ~ x, data = dd, prior = p1))
  at("frmtmb frm(p3 class=foo)", frm(y ~ x, data = dd, prior = p3))
  at("frmtmb frm(p2)", frm(y ~ x, data = dd, prior = p2))
}
if (arm == "lane") {
  suppressMessages(library(brms))
  at("brms set_prior(class=sigma)",
     b1 <- brms::set_prior("dirichlet(2)", class = "sigma"))
  at("brms validate_prior(b1)",
     brms::validate_prior(b1, y ~ x, data = dd))
  at("brms stancode(b1)", brms::stancode(y ~ x, data = dd, prior = b1))
  at("brms set_prior(class=b)",
     b2 <- brms::set_prior("dirichlet(2)", class = "b"))
  at("brms stancode(b2)", brms::stancode(y ~ x, data = dd, prior = b2))
}
