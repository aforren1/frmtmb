# Exploration only: what frmtmb's frame, family and prior objects hold,
# so the shims in tests/testthat/helper-brms-suite.R map brms's names
# onto something that exists. Run: Rscript dev/brmsport-probe.R
lib <- Sys.getenv("BRMSPORT_LIB", "C:/Users/adf44/source/r/brmsport-lib")
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "\n")
ex <- getNamespaceExports("frmtmb")
cat(length(ex), "exports\n")
want <- c("brmsfamily", "default_prior", "validate_prior", "as.brmsprior",
          "empty_prior", "update_adterms", "brmsterms", "make_conditions",
          "prior_string", "prior_", "mvbind", "set_rescor", "lf", "nlf",
          "xbeta", "hurdle_negbinomial", "is.brmsformula", "ngrps",
          "conditional_effects", "conditional_smooths", "pp_check",
          "hypothesis", "variables", "ndraws", "nchains", "niterations",
          "posterior_summary", "model_weights", "nsamples", "parnames",
          "posterior_samples", "rename_pars", "standata", "stancode",
          "prior_summary", "autocor", "theme_black", "probit_approx",
          "inv_link", "log_lik", "residuals", "mixture", "thres", "me",
          "gr", "mm", "mmc", "cs", "mo", "arma", "ar", "ma", "gp", "s",
          "t2", "rate", "subset", "index", "mi", "cens", "trials", "se",
          "weights", "cat", "dec", "vreal", "posterior_epred",
          "posterior_linpred", "posterior_predict", "predictive_error",
          "pp_mixture", "as_draws", "prior_draws", "logm1", "expm1")
print(setNames(want %in% ex, want))

set.seed(1)
d <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))
fr <- frm(y ~ x, d, dry_run = "frame")
cat("\nframe class:", class(fr), "\n")
str(fr, max.level = 1, give.attr = FALSE)
fr2 <- frm(y | cens(c1) ~ 1,
           data.frame(y = rnorm(9), c1 = rep(-1:1, 3)), dry_run = "frame")
str(fr2, max.level = 2, give.attr = FALSE)
cat("\nfamily fields:\n")
print(names(student()))
print(student()$link)
print(class(set_prior("normal(0, 2)", class = c("b", "sd"))))
