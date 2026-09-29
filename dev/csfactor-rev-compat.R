# Claims 3 and 5: the frm_sample() narrowing, and the stored-layout
# consumers (frm_simulate(newparams=), update(start=), frm_compat(),
# frmtmb.sample on a cs() factor fit, variables() on a draws object).
#   Rscript dev/csfactor-rev-compat.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("== build", TAG, ": frmtmb", as.character(packageVersion("frmtmb")),
    "frmtmb.sample", as.character(packageVersion("frmtmb.sample")),
    "from", dirname(system.file(package = "frmtmb.sample")), "\n")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                           1L, 300L)
go <- function(lab, expr) {
  cat("\n-- ", lab, " --\n", sep = "")
  print(tryCatch(expr, error = function(e) paste("ERROR:", short(e))))
}

# the worker's construction, verbatim: seed 405, n = 250
set.seed(405)
n <- 250
fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(fc)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n)
d <- data.frame(fc = fc,
                yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
d$x <- rnorm(n)

fit <- frm(bf(yo ~ cs(fc)), family = sratio(), data = d)
cat("ML fit logLik = ", sprintf("%.7f", as.numeric(logLik(fit))),
    "  npar = ", length(fit$opt$par), "\n", sep = "")
cat("par_template components: ",
    paste(names(fit$estimates), collapse = ", "), "\n", sep = "")
for (nmc in grep("^bcs", names(fit$estimates), value = TRUE)) {
  cat("  ", nmc, " = ", paste(sprintf("%.6f", fit$estimates[[nmc]]),
                              collapse = " "), "\n", sep = "")
}

cat("\n########## claim 3: frm_sample() on y ~ x + cs(x) ##########\n")
go("frm_sample(yo ~ x + cs(x))",
   frm_sample(bf(yo ~ x + cs(x)), family = sratio(), data = d,
              chains = 1L, iter = 200L, refresh = 0L, seed = 3L))
go("frm(yo ~ x + cs(x)) for the same message",
   frm(bf(yo ~ x + cs(x)), family = sratio(), data = d))
go("frm_sample with a proper prior on b and bcs",
   frm_sample(bf(yo ~ x + cs(x)), family = sratio(), data = d,
              prior = set_prior("normal(0, 1)", class = "b"),
              chains = 1L, iter = 200L, refresh = 0L, seed = 3L))

cat("\n########## claim 5: stored layouts ##########\n")
go("variables(fit)", variables(fit))
go("frm_compat() rows mentioning cs", {
  tb <- frm_compat()
  tb[grepl("cs", tb[[1L]]) | grepl("cs", tb[[2L]]), , drop = FALSE]
})
go("frm_simulate(newparams=) with the fit's own cs components", {
  np <- fit$estimates[c("b", "tau_raw",
                        grep("^bcs", names(fit$estimates), value = TRUE))]
  np <- np[lengths(np) > 0L]
  s <- frm_simulate(bf(yo ~ cs(fc)), family = sratio(), data = d,
                    newparams = np, nsim = 2, seed = 4L)
  c(n = length(s), first_table = paste(table(s[[1L]]), collapse = "/"))
})
go("frm_simulate(newparams=) with a WRONG cs length", {
  np <- fit$estimates[c("b", "tau_raw",
                        grep("^bcs", names(fit$estimates), value = TRUE))]
  np <- np[lengths(np) > 0L]
  k <- grep("^bcs", names(np), value = TRUE)[1L]
  np[[k]] <- c(np[[k]], 0)
  frm_simulate(bf(yo ~ cs(fc)), family = sratio(), data = d,
               newparams = np, nsim = 1, seed = 4L)
})
go("update(start = the fit's estimates)", {
  u <- update(fit, start = fit$estimates)
  c(logLik = sprintf("%.7f", as.numeric(logLik(u))),
    npar = length(u$opt$par))
})
go("update(formula. = ~ x + cs(fc))",
   update(fit, formula. = yo ~ x + cs(fc)))
go("update(formula. = ~ fc + cs(fc)) must be refused",
   update(fit, formula. = yo ~ fc + cs(fc)))

cat("\n########## frmtmb.sample on the cs() factor fit ##########\n")
ds <- tryCatch(frm_sample(bf(yo ~ cs(fc)), family = sratio(), data = d,
                          chains = 2L, iter = 600L, refresh = 0L,
                          seed = 7L),
               error = function(e) structure(list(m = short(e)),
                                             class = "revfail"))
if (inherits(ds, "revfail")) {
  cat("frm_sample FAILED: ", ds$m, "\n")
} else {
  go("variables(ds) bcs rows", grep("^bcs", variables(ds), value = TRUE))
  go("variables(ds) all", variables(ds))
  go("fixef(ds) rows", rownames(fixef(ds)))
  go("posterior_epred at factor('c'), posterior mean",
     round(apply(posterior_epred(ds, newdata = data.frame(
       fc = factor("c"))), c(2L, 3L), mean), 4))
  go("posterior_epred over the three levels",
     round(apply(posterior_epred(ds, newdata = data.frame(
       fc = factor(c("a", "b", "c")))), c(2L, 3L), mean), 4))
  go("ML fit at the same rows",
     round(frm_linpred(fit, newdata = data.frame(
       fc = factor(c("a", "b", "c"))), type = "response"), 4))
}
cat("\nDONE ", TAG, "\n")
