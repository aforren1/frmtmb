# Release review of 0.74.0, item 6: the claims of core's NEWS, probed on
# rellib-r12 (frmtmb 0.74.0) and, where a bullet says "before" or
# "through 0.73.0", on rellib-r11.
#   Rscript dev/review-rel074-news.R r12|r11
a <- commandArgs(trailingOnly = TRUE)
lib <- if (a[1] == "r11") "C:/Users/adf44/source/r/rellib-r11" else
  "C:/Users/adf44/source/r/rellib-r12"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "\n")
p <- function(label, expr) {
  r <- tryCatch({
    v <- withCallingHandlers(expr, warning = function(w) {
      cat("   warning:", substr(conditionMessage(w), 1, 140), "\n")
      invokeRestart("muffleWarning")})
    if (inherits(v, "frmtmb_fit")) sprintf("FIT logLik %.8f", as.numeric(logLik(v)))
    else paste("OK", paste(utils::capture.output(str(v, max.level = 0, give.attr = FALSE)), collapse = " "))
  }, error = function(e) paste("ERROR", substr(gsub("\\s+", " ", conditionMessage(e)), 1, 260)))
  cat(sprintf("%-44s %s\n", label, r))
}
set.seed(11)
d <- data.frame(g = factor(rep(1:12, each = 8)), x = rnorm(96))
d$y <- 1 + 0.5 * d$x + rnorm(12, 0, 0.8)[d$g] + rnorm(96, 0, 0.5)
foo <- function(x) x
p("(1 | g)", frm(bf(y ~ x + (1 | g)), data = d))
p("(1 | gr(g))", frm(bf(y ~ x + (1 | gr(g))), data = d))
p("foo(1 | g)", frm(bf(y ~ x + foo(1 | g)), data = d))
p("us_t(1 | g)", frm(bf(y ~ x + us_t(1 | g)), data = d))
p("set_prior(dirichlet, class = sigma)",
  set_prior("dirichlet(1)", class = "sigma"))
p("frm(prior = dirichlet on sigma)",
  frm(bf(y ~ x), data = d, prior = set_prior("dirichlet(1)", class = "sigma")))
p("validate_prior(dirichlet on sigma)",
  validate_prior(set_prior("dirichlet(1)", class = "sigma"), bf(y ~ x),
                 data = d))
p("set_prior(dirichlet, class = b)", set_prior("dirichlet(1)", class = "b"))
if (a[1] == "r12") {
  p("us_sigma(theta length 2, d = 2)", us_sigma(c(0, 0), 2))
  p("us_sigma(theta length 3, d = 2)", us_sigma(c(0, 0, 0.3), 2))
  p("us_theta(us_sigma(th, 3)) == th",
    { th <- c(0.1, -0.2, 0.3, 0.2, -0.1, 0.4); all.equal(us_theta(us_sigma(th, 3)), th) })
  dummy <- function(...) stop("unused")
  for (nm in c("gaussian", "dec", "poly", "ns", "log", "car")) {
    p(paste0("register_covstruct(\"", nm, "\")"),
      frmtmb_register_covstruct(nm, kind = "bar", n_theta = function(...) 1L,
                                lpdf = dummy))
  }
  print(args(frmtmb_register_covstruct))
}
