## Reviewer: read brms's own bodies, independently of the lane's script.
LIB <- "C:/Users/adf44/source/r/bmrev-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
cat("posterior", format(packageVersion("posterior")), "\n\n")

show1 <- function(nm) {
  f <- tryCatch(get(nm, envir = asNamespace("brms")),
                error = function(e) NULL)
  cat("========== brms:::", nm, " ==========\n", sep = "")
  if (is.null(f)) { cat("  ABSENT\n\n"); return(invisible()) }
  cat("formals: ", paste(names(formals(f)), collapse = " "), "\n")
  print(f)
  cat("\n")
}
for (nm in c("rhat.brmsfit", "neff_ratio.brmsfit",
             "nuts_params.brmsfit", "log_posterior.brmsfit",
             "extract_pars", "as.mcmc.brmsfit",
             "posterior_interval.brmsfit")) show1(nm)

cat("========== extract_pars refusal ==========\n")
ap <- c("b_Intercept", "b_x", "sigma", "lp__")
for (lbl in list(list("TRUE", TRUE), list("0.9", 0.9),
                 list("NA", NA), list('"^b"', "^b"),
                 list("NULL", NULL), list("1L", 1L))) {
  v <- tryCatch(brms:::extract_pars(lbl[[2]], all_pars = ap,
                                    fixed = FALSE),
                error = function(e) paste0("ERROR: ",
                                           conditionMessage(e)))
  cat(sprintf("%-8s -> %s\n", lbl[[1]],
              paste(utils::capture.output(str(v)), collapse = " ")))
}
cat("DONE\n")
