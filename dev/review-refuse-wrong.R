# Item 1: does a deliberately wrong model still FAIL through the edited
# translator? Each mutant is applied in this process, one at a time, to the
# namespace or to the translator's own text. Seed 9102.
source("dev/review-refuse-lib.R")
set.seed(9102)
n <- 150
x <- rnorm(n)
ds <- list()
ds$Gamma <- data.frame(x = x, y = rgamma(n, 2.5, 2.5 / exp(0.3 + 0.5 * x)))
ds$Beta <- data.frame(x = x, y = rbeta(n, plogis(0.2 + 0.6 * x) * 6,
                                      (1 - plogis(0.2 + 0.6 * x)) * 6))
ds$von_mises <- data.frame(x = x, y = ((2 * atan(0.3 + 0.5 * x) +
                           rnorm(n, 0, 0.5) + pi) %% (2 * pi)) - pi)
for (k in names(ds)) ds[[k]]$y[seq(7, 147, by = 14)] <- NA
fam_txt <- c(Gamma = "Gamma(link = 'log')", Beta = "Beta()",
             von_mises = "von_mises()")
one <- function(fam, label, extra = list()) {
  sp <- c(list(bf = "bf(y | mi() ~ x)", family = fam_txt[[fam]],
               data = ds[[fam]]), extra)
  r <- tr(meas(sp))
  if (is_err(r)) {
    pr("%-34s %-10s ERROR %s", label, fam, msg(r))
    return(invisible())
  }
  pr("%-34s %-10s verdict=%-14s const=%s pts=[%s] grad=%s err=%s", label, fam,
     r$verdict, format(r$l1_const %||% NA, digits = 3),
     paste(format(unlist(r$l1_extra_const) %||% NA, digits = 3), collapse = ","),
     format(r$l1_grad %||% NA, digits = 3),
     substr(r$l1_error %||% r$frm_error %||% "", 1, 80))
}
ns <- asNamespace("frmtmb")
orig_j <- ns$mi_chart_logjac
orig_i <- ns$mi_chart_inv
orig_m <- ns$mi_chart_maps
for (fam in names(ds)) one(fam, "control: lane's code")

# M1: the objective drops the Jacobian, the translator still removes it
mutj <- function(fn) assignInNamespace("mi_chart_logjac", fn, "frmtmb")
assignInNamespace("mi_chart_maps", function(chart) FALSE, "frmtmb")
for (fam in names(ds)) one(fam, "M1 objective has no Jacobian")
assignInNamespace("mi_chart_maps", orig_m, "frmtmb")

# M2: a wrong Jacobian that the translator reads too: the constant cancels by
# construction, so only the inner gradient can see it
mutj(function(u, chart) 2 * orig_j(u, chart))
for (fam in names(ds)) one(fam, "M2 Jacobian doubled on both sides")
mutj(orig_j)

# M3: a Jacobian off by a parameter-free constant, which no L1 sees and which
# changes nothing about the model
mutj(function(u, chart) orig_j(u, chart) + 1)
for (fam in names(ds)) one(fam, "M3 Jacobian + 1 per value")
mutj(orig_j)

# T1, T2: the translator alone. Its text is rewritten; the objective stands
h <- readLines("tests/testthat/helper-brms.R")
hj <- h[!grepl("jac <- jac - sum(frmtmb:::mi_chart_logjac", h, fixed = TRUE)]
stopifnot(length(hj) == length(h) - 1L)
tf <- tempfile(fileext = ".R")
writeLines(hj, tf)
sys.source(tf, envir = E)
for (fam in names(ds)) one(fam, "T1 translator keeps the Jacobian")
hv <- h
at <- grep("v <- frmtmb:::mi_chart_inv(u, mm$chart)", h, fixed = TRUE)
stopifnot(length(at) == 1L)
hv[at] <- "      v <- u"
writeLines(hv, tf)
sys.source(tf, envir = E)
for (fam in names(ds)) one(fam, "T2 translator passes u as the value")
tf2 <- "tests/testthat/helper-brms.R"
sys.source(tf2, envir = E)

# W1: another family on the brms side: the same frm() fit, a different program
one("Gamma", "W1 brms is weibull, frm is Gamma",
    list(family_brms = "weibull()", family_frm = "Gamma(link = 'log')"))
one("Beta", "W2 brms probit, frm logit",
    list(family_brms = "Beta(link = 'probit')", family_frm = "Beta()"))
