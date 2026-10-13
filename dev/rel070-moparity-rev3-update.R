# Reviewer, final check (B2): slot matching of a prior table's rows
# against stored specifications spelled differently (class sigma, dpar
# rows, a sigma Intercept row against mu's), and a stored legacy list
# plus a table.
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(11)
dd <- data.frame(x = rnorm(80), z = rnorm(80))
dd$y <- rnorm(80, 1 + 0.5 * dd$x, exp(0.2 + 0.3 * dd$z))
rows <- function(f) {
  if (is.null(f$prior)) return("<none>")
  if (!inherits(f$prior, "frmtmb_priorlist")) return(paste("named list:", paste(names(f$prior), collapse = ",")))
  d <- as.data.frame(f$prior)
  paste(paste0(d$class, "/", d$coef, "/", d$dpar, ":", d$prior), collapse = "; ")
}
try1 <- function(lab, expr) {
  r <- tryCatch(suppressMessages(expr), error = function(e) e)
  cat("==", lab, ":", if (inherits(r, "error")) paste("ERROR", substr(conditionMessage(r), 1, 150)) else rows(r), "\n")
}
fo <- bf(y ~ x, sigma ~ z)
p0 <- set_prior("normal(0, 0.5)", class = "b", coef = "z", dpar = "sigma") +
  set_prior("normal(0, 10)", class = "Intercept") +
  set_prior("normal(0, 2)", class = "Intercept", dpar = "sigma")
f0 <- frm(fo, family = gaussian(), data = dd, prior = p0)
cat("stored:", rows(f0), "\n")
tb <- brms::default_prior(brms::bf(y ~ x, sigma ~ z), data = dd)
cat("brms table rows:", paste(paste0(tb$class, "/", tb$coef, "/", tb$dpar), collapse = "; "), "\n")
one <- function(cls, coef = "", dpar = "") {
  r <- tb[tb$class == cls & tb$coef == coef & tb$dpar == dpar, ]
  r$prior <- ""
  r
}
try1("empty row b z dpar sigma", update(f0, prior = one("b", "z", "sigma")))
try1("empty row Intercept dpar sigma", update(f0, prior = one("Intercept", "", "sigma")))
try1("empty row Intercept (mu)", update(f0, prior = one("Intercept")))
try1("empty row b (mu, class-wide)", update(f0, prior = one("b")))
# a stored class = "sigma" spelling against a table row class sigma
dd2 <- dd
f1 <- frm(bf(y ~ x), family = gaussian(), data = dd2,
          prior = set_prior("exponential(3)", class = "sigma") +
            set_prior("normal(0, 10)", class = "Intercept"))
cat("stored f1:", rows(f1), "\n")
tb1 <- brms::default_prior(y ~ x, data = dd2)
s1 <- tb1[tb1$class == "sigma", ]; s1$prior <- ""
try1("f1, empty sigma row", update(f1, prior = s1))
s2 <- tb1[tb1$class == "Intercept", ]; s2$prior <- ""
try1("f1, empty Intercept row", update(f1, prior = s2))
# legacy named list stored, then a table
f2 <- frm(bf(y ~ x), family = gaussian(), data = dd2,
          prior = list(beta = prior_normal(0, 0.05)))
try1("legacy list stored + table", update(f2, prior = s1))
try1("legacy list stored + zero-row table", update(f2, prior = brms::empty_prior()))
