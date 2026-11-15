# Re-check of lane priors2 after fix round 2, item 3: a fit beside an
# open end whose box run raises "NA/NaN function evaluation" and whose box
# result is returned; item 4: cauchit bounds outside (0, 1).
# Rscript priors2-short-r2b.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
show <- function(lab, expr) {
  ws <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, paste0(class(w)[1], ": ", substr(conditionMessage(w), 1, 80)))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    cat(lab, ": ERROR ", substr(conditionMessage(r), 1, 150), "\n", sep = "")
  } else {
    cat(lab, ": code ", r$opt$convergence, " hex ",
        paste(sprintf("%a", c(r$opt$objective, r$opt$par)), collapse = ","),
        "\n", sep = "")
  }
  for (w in ws) cat("    [warning]", w, "\n")
  invisible(r)
}
sp <- set_prior
set.seed(11)
pd <- data.frame(x = runif(80, 0.5, 2), w = runif(80, 0, 3))
pd$y <- rpois(80, 1.5 * pd$x)
pd$y2 <- rpois(80, 0.2 + 1.5 * pd$x + 0.5 * pd$w)
pd$ye <- rexp(80, 1 / (0.5 + 0.8 * pd$x))
pd$one <- 1
for (pr in c("lognormal(0, 1)", "gamma(2, 2)", "inv_gamma(3, 1)",
             "beta(2, 2)")) {
  ub <- if (startsWith(pr, "beta")) 3 else NA
  show(sprintf("poisson(identity) y2 ~ x + w, %s on w lb=0", pr),
       frm(y2 ~ x + w, data = pd, family = poisson(link = "identity"),
           prior = sp(pr, class = "b", coef = "w", lb = 0, ub = ub)))
  show(sprintf("poisson(identity) y2 ~ x + w, %s on x lb=0", pr),
       frm(y2 ~ x + w, data = pd, family = poisson(link = "identity"),
           prior = sp(pr, class = "b", coef = "x", lb = 0, ub = ub)))
  show(sprintf("exponential(identity) ye ~ x, %s on x lb=0", pr),
       frm(ye ~ x, data = pd, family = exponential(link = "identity"),
           prior = sp(pr, class = "b", coef = "x", lb = 0, ub = ub)))
  show(sprintf("exponential(identity) ye ~ 0 + one, %s lb=0", pr),
       frm(ye ~ 0 + one, data = pd, family = exponential(link = "identity"),
           prior = sp(pr, class = "b", coef = "one", lb = 0, ub = ub)))
  show(sprintf("poisson(identity) y ~ 0 + x, %s lb=0, start 1e-3", pr),
       frm(y ~ 0 + x, data = pd, family = poisson(link = "identity"),
           prior = sp(pr, class = "b", coef = "x", lb = 0, ub = ub),
           start = list(beta = 1e-3)))
}
cat("\n## cauchit: a bound outside (0, 1) through tan()\n")
pm <- list(scale = "link", link = frmtmb:::get_link("cauchit"))
lk <- frmtmb:::get_link("cauchit")
for (s in c("lb", "ub")) for (v in c(-0.5, -1, 1.5, 2, 2.5)) {
  r <- tryCatch(frmtmb:::internal_bound(v, pm, s), error = function(e) NA)
  cat(sprintf("  cauchit %s = %4s -> b = %s, linkinv(b) = %s\n", s, v,
              format(r, digits = 6),
              if (is.na(r)) "refused" else format(lk$linkinv(r), digits = 6)))
}
set.seed(1)
dz <- data.frame(z = rnorm(200))
dz$cnt <- rpois(200, 2) * rbinom(200, 1, 0.8)
fam <- tryCatch(zero_inflated_poisson(link_zi = "cauchit"),
                error = function(e) NULL)
if (!is.null(fam)) {
  show("zip cauchit zi, no prior", frm(cnt ~ z, data = dz, family = fam))
  show("zip cauchit zi, lb = -0.5",
       frm(cnt ~ z, data = dz, family = fam,
           prior = sp("", class = "zi", lb = -0.5)))
  show("zip cauchit zi, ub = 1.5",
       frm(cnt ~ z, data = dz, family = fam,
           prior = sp("", class = "zi", ub = 1.5)))
} else cat("zero_inflated_poisson(link_zi = 'cauchit') not available\n")
