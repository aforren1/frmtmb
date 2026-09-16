## Reviewer: the last few probes.
LIB <- "C:/Users/adf44/source/r/bmrev-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample)); q(library(posterior))
ds <- readRDS("dev/stan-cache/bmrev-draws.rds")
go <- function(lbl, e) {
  v <- try(q(e), silent = TRUE)
  if (inherits(v, "try-error"))
    cat(sprintf("%-46s ERROR: %s\n", lbl,
                sub("\n.*$", "",
                    conditionMessage(attr(v, "condition")))))
  else cat(sprintf("%-46s %s [%s]\n", lbl,
                   paste(class(v), collapse = "/"),
                   paste(dim(v) %||% length(v), collapse = "x")))
  invisible(v)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

cat("== nuts_params / log_posterior really delegate ==\n")
np <- q(nuts_params(ds))
cat("nuts_params(ds) class/dim: ", paste(class(np), collapse = "/"), " ",
    paste(dim(np), collapse = "x"), "\n")
cat("Parameter levels: ",
    paste(levels(factor(np$Parameter)), collapse = " "), "\n")
bn <- q(bayesplot::nuts_params(ds$stanfit))
cat("identical to bayesplot::nuts_params(ds$stanfit): ",
    identical(np, bn), "\n")
np2 <- q(nuts_params(ds, "stepsize__"))
cat("nuts_params(ds, \"stepsize__\") positional -> levels: ",
    paste(levels(factor(np2$Parameter)), collapse = " "), "\n")
lp <- q(log_posterior(ds))
cat("log_posterior(ds) identical to bayesplot's on stanfit: ",
    identical(lp, q(bayesplot::log_posterior(ds$stanfit))), "\n")
cat("dim: ", paste(dim(lp), collapse = "x"), "\n")

cat("\n== bayes_R2 slot 4: brms's `robust`, ours `probs` ==\n")
b0 <- q(bayes_R2(ds))
cat("bayes_R2(ds) colnames:                 ",
    paste(colnames(b0), collapse = " "), "\n")
b1 <- q(bayes_R2(ds, NULL, TRUE, TRUE))
cat("bayes_R2(ds, NULL, TRUE, TRUE) cols:   ",
    paste(colnames(b1), collapse = " "), "\n")
cat("  -> brms reads slot 4 as robust=TRUE (median/MAD);",
    "here it is probs=TRUE\n")
cat("  answers silently rather than erroring: ",
    !inherits(try(q(bayes_R2(ds, NULL, TRUE, TRUE)), silent = TRUE),
              "try-error"), "\n")

cat("\n== pairs() error quality ==\n")
go("pairs(ds, \"nosuchvariable\")", pairs(ds, "nosuchvariable"))
go("pairs(ds, \"^x$\")", pairs(ds, "^x$"))
go("mcmc_plot(ds, variable = \"nosuchvariable\")",
   mcmc_plot(ds, variable = "nosuchvariable"))

cat("\n== summary(ds, prob) slot: brms's position 3 ==\n")
s0 <- q(summary(ds)); s1 <- q(summary(ds, NULL, 0.5))
cat("summary(ds) cols:           ", paste(colnames(s0), collapse = " "),
    "\n")
cat("summary(ds, NULL, 0.5) cols:", paste(colnames(s1), collapse = " "),
    "\n")
cat("identical: ", identical(s0, s1), "\n")
cat("brms: summary(fit, priors, prob = 0.5) gives l-50% / u-50% CI\n")
cat("brms summary column names are Rhat, Bulk_ESS, Tail_ESS\n")
cat("DONE\n")
