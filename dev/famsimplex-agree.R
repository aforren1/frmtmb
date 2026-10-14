# ML fits of the simplex families against brms 2.23.0 fits of the same
# data (dev/famsimplex-brms-fits2.R, 4 chains x 1000 post-warmup draws,
# seed 1): (posterior median - ML estimate) / posterior sd, per brms
# variable. ML values are read under brms's names (hyp_env_vals()).
# Output: dev/famsimplex-agree.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex"
forms <- list(dir = bf(Y ~ x, phi ~ z), dirre = bf(Y ~ x + (1 | ID | g)),
              ln = bf(Y ~ x, sigmab ~ z), lnre = bf(Y ~ x + (1 | ID | g)),
              dir2 = bf(Y ~ x))
fams <- list(dir = dirichlet(), dirre = dirichlet(),
             ln = logistic_normal(), lnre = logistic_normal(),
             dir2 = brmsfamily("dirichlet2"))
all <- list()
for (nm in names(fams)) {
  b <- readRDS(file.path(wt, "dev", paste0("famsimplex-brms-", nm, ".rds")))
  fit <- frm(forms[[nm]], family = fams[[nm]], data = b$data)
  vo <- ns$hyp_vals_only(fit)
  ml <- unlist(ns$hyp_env_vals(fit, vo$vals, vo$comp))
  s <- b$summ
  s <- s[!grepl("^Intercept_", s$variable), ]
  miss <- setdiff(s$variable, names(ml))
  if (length(miss)) cat(nm, "brms variables with no ML value:", miss, "\n")
  s <- s[s$variable %in% names(ml), ]
  d <- data.frame(model = nm, variable = s$variable, brms_median = s$median,
                  brms_sd = s$sd, ml = unname(ml[s$variable]),
                  rhat = s$rhat)
  d$z <- (d$brms_median - d$ml) / d$brms_sd
  all[[nm]] <- d
  cat("\n==", nm, ": converged", fit$opt$convergence == 0, " logLik",
      format(as.numeric(logLik(fit)), digits = 10), "\n")
  print(d, digits = 4, row.names = FALSE)
}
out <- do.call(rbind, all)
cat("\nSUMMARY |z| per model\n")
print(aggregate(abs(z) ~ model, out, function(v) c(max = max(v),
                                                    mean = mean(v))),
      digits = 3)
cat("variables compared:", nrow(out), " max |z|:",
    format(max(abs(out$z)), digits = 3), "\n")
