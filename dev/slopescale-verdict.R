# Lane slopescale, item 1018, brief step 3: the boundary verdict on a
# short stop. With autoscale = FALSE a random slope on a column of
# spread 1e-4 or 1e-5 still stalls at its start; what does the fit say?
# The triage saw "sd_g__x is at zero" at 1e-5, where the sd is 1.0 and
# the maximum is near 5e4. Counts every condition class and whether a
# boundary message names the slope's sd.
# Usage: Rscript slopescale-verdict.R before|after
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")),
    find.package("frmtmb"), "arm", arm, "\n")
rows <- list()
first <- character(0)
for (fam in c("gau", "pois")) for (form in c("0", "1"))
  for (s in c(1e-4, 1e-5)) for (seed in 1:6) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x0 <- rnorm(300)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * x0
  y <- if (fam == "gau") eta + rnorm(300) else rpois(300, exp(eta))
  d <- data.frame(g = g, x = x0 * s, y = y)
  fo <- if (form == "0") y ~ 1 + (0 + x | g) else y ~ 1 + (1 + x | g)
  fm <- if (fam == "gau") gaussian() else poisson()
  cl <- character(0)
  txt <- character(0)
  f <- withCallingHandlers(
    frm(fo, family = fm, data = d,
        control = frmtmb_control(autoscale = FALSE)),
    warning = function(w) {
      cl <<- c(cl, paste0("W:", class(w)[1]))
      txt <<- c(txt, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cl <<- c(cl, paste0("M:", class(m)[1]))
      txt <<- c(txt, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  zero <- any(grepl("sd_g__x is at zero", txt, fixed = TRUE))
  key <- paste(fam, form, s)
  if (!key %in% first) {
    first <- c(first, key)
    cat("\n==", key, "seed", seed, "\n")
    cat(paste0("  ", substr(txt, 1, 300)), sep = "\n")
  }
  rows[[length(rows) + 1L]] <- data.frame(
    fam = fam, form = form, s = s, seed = seed,
    sd_x = unname(exp(f$opt$par[grep("theta", names(f$opt$par))][
      if (form == "0") 1 else 2])),
    conv = f$opt$convergence, n_warn = sum(startsWith(cl, "W:")),
    at_zero = zero, classes = paste(unique(cl), collapse = " "))
}
r <- do.call(rbind, rows)
options(width = 200)
print(r, digits = 4)
cat("\nfits:", nrow(r), "; 'sd_g__x is at zero':", sum(r$at_zero),
    "; with no warning:", sum(r$n_warn == 0), "\n")
