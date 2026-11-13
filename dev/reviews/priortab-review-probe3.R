# Reviewer of lane priortab: 4266 heuristic, interior MAP under a closed
# support box, frm_sample(fit) after a bounds-only prior.
# PRIORTAB_LIB unset = lane library; PRIORTAB_LIB=none = rellib-r20.
lib <- Sys.getenv("PRIORTAB_LIB", "C:/Users/adf44/source/r/wt-priortab-lib")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
options(width = 150, digits = 10)
show <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", substr(gsub("\n", " ", conditionMessage(w)), 1, 140),
               "\n")
           invokeRestart("muffleWarning")
         },
         message = function(m) invokeRestart("muffleMessage")),
       error = function(e) {
         cat("ERROR: ", substr(conditionMessage(e), 1, 300), "\n")
         NULL
       })
  if (!is.null(r)) print(r)
  invisible(r)
}
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$yneg <- d$y - 0.8 * d$z
fml <- frm(y ~ z + (1 | g), data = d)
eff <- function(pl) {
  rp <- frmtmb.sample:::sample_resolve_priors(fml, pl, base = NULL,
                                              defaults = TRUE)
  x <- as.data.frame(rp$effective)
  x[, intersect(c("prior", "class", "coef", "group", "lb"), names(x))]
}
show("H1 coerced fit table, sigma row set to brms's default", {
  tab <- as.data.frame(default_prior(y ~ z + (1 | g), data = d))
  tab$prior[tab$class == "sigma"] <- "student_t(3, 0, 2.5)"
  eff(as.brmsprior(tab))
})
show("H2 coerced fit table, untouched", {
  eff(as.brmsprior(as.data.frame(default_prior(y ~ z + (1 | g), data = d))))
})
show("H3 coerced sample table, untouched", {
  eff(as.brmsprior(as.data.frame(default_prior(y ~ z + (1 | g), data = d,
                                               route = "sample"))))
})
show("H4 coerced sample table, Intercept edited to (flat)", {
  tab <- as.data.frame(default_prior(y ~ z + (1 | g), data = d,
                                     route = "sample"))
  tab$prior[tab$class == "Intercept"] <- "(flat)"
  eff(as.brmsprior(tab))
})

for (dn in c("lognormal(0, 1)", "beta(2, 2)", "gamma(2, 2)")) {
  show(paste("M", dn, "no bound; yneg"), {
    f <- frm(yneg ~ z, data = d, prior = set_prior(dn, class = "b",
                                                  coef = "z"))
    c(conv = f$opt$convergence, obj = f$opt$objective,
      z = unname(fixef(f)["z", "Estimate"]))
  })
  show(paste("M", dn, "lb = 1e-8; yneg"), {
    f <- frm(yneg ~ z, data = d, prior = set_prior(dn, class = "b",
                                                  coef = "z", lb = 1e-8,
                                                  ub = if (grepl("beta", dn)) 1))
    c(conv = f$opt$convergence, obj = f$opt$objective,
      z = unname(fixef(f)["z", "Estimate"]))
  })
}
cat("\ndone\n")
