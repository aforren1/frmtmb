## What brms does with the two positional calls that answer silently
## here: as.mcmc(x, TRUE) is pars = TRUE and posterior_interval(x, 0.9)
## is pars = 0.9. Run: Rscript dev/brmsmatch-pars.R
.libPaths(c("C:/Users/adf44/source/r/brmsmatch-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
loadNamespace("brms")

print(brms:::extract_pars)
print(brms:::posterior_interval.brmsfit)

all_pars <- c("Intercept", "x", "sigma", "b[1]", "lp__")
for (p in list(NA, TRUE, 0.9, "x", c("Intercept", "x"))) {
  cat("\npars = ", paste(deparse(p), collapse = ""), "\n", sep = "")
  r <- try(brms:::extract_pars(p, all_pars = all_pars, fixed = FALSE),
           silent = TRUE)
  print(r)
}

## ROUND 2. brms spells TWO different selectors `pars`, and the round-1
## review found this lane had given rhat()/neff_ratio() the wrong one.
## brms:::rhat.brmsfit passes `variable = pars` to as_draws_array()
## rather than through extract_pars(), so NULL is its default and a
## string is an EXACT name. The five probed values, here:
.libPaths(c("C:/Users/adf44/source/r/brmsmatch-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
qq <- function(e) suppressWarnings(suppressMessages(e))
qq(library(frmtmb)); qq(library(frmtmb.sample))
ds <- readRDS("dev/stan-cache/brmsmatch-draws.rds")
cat("\n== rhat(ds, <pars>) after the round-2 fix ==\n")
for (v in list(quote(rhat(ds)), quote(rhat(ds, "x")),
               quote(rhat(ds, "^x$")), quote(rhat(ds, NULL)),
               quote(rhat(ds, NA)),
               quote(rhat(ds, "^x$", regex = TRUE)))) {
  r <- try(qq(eval(v)), silent = TRUE)
  cat(sprintf("%-34s %s\n", paste(deparse(v), collapse = ""),
              if (inherits(r, "try-error")) {
                paste0("ERROR: ", substr(sub("\n.*$", "",
                  sub("^Error[^:]*: ", "", as.character(r))), 1, 44))
              } else {
                paste0(length(r), " variable(s): ",
                       paste(utils::head(names(r), 3), collapse = " "))
              }))
}
cat("DONE\n")
