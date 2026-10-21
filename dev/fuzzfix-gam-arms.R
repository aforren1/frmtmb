# Minimize the Gamma + weights + ar1 + mo_int + shape ~ z, profile spec
# (fuzz seed 20502220) one component at a time.
#   Rscript dev/fuzzfix-gam-arms.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
d <- readRDS("dev/fuzzfix-out/data-gam.rds")
G <- Gamma(link = "log")
full_rhs <- "1 + x + mo(xo) * z + ar1(tt + 0 | g)"
arm <- function(lhs = "y | weights(w)", rhs = full_rhs, shape = TRUE,
                profile = TRUE) {
  list(lhs = lhs, rhs = rhs, shape = shape, profile = profile)
}
arms <- list(
  full = arm(),
  no_weights = arm(lhs = "y"),
  no_profile = arm(profile = FALSE),
  no_shape = arm(shape = FALSE),
  re_ri = arm(rhs = "1 + x + mo(xo) * z + (1 | g)"),
  re_none = arm(rhs = "1 + x + mo(xo) * z"),
  mo_main = arm(rhs = "1 + x + mo(xo) + z + ar1(tt + 0 | g)"),
  mo_none = arm(rhs = "1 + x + xo * z + ar1(tt + 0 | g)"),
  min_pro = arm(lhs = "y", rhs = "1 + x + mo(xo) * z", shape = FALSE),
  min_pro_mo = arm(lhs = "y", rhs = "1 + mo(xo)", shape = FALSE),
  min_pro_nomo = arm(lhs = "y", rhs = "1 + x + xo * z", shape = FALSE)
)
for (nm in names(arms)) {
  ar <- arms[[nm]]
  fo <- stats::as.formula(paste(ar$lhs, "~", ar$rhs))
  bfo <- if (ar$shape) bf(fo, shape ~ 1 + z) else bf(fo)
  ctl <- frmtmb_control(profile = ar$profile)
  w <- character()
  f <- withCallingHandlers(
    tryCatch(frm(bfo, family = G, data = d, control = ctl),
             error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  cat(sprintf("ARM %-12s %s | warnings: %d %s\n", nm,
              if (inherits(f, "error")) paste("ERROR",
                substr(conditionMessage(f), 1, 110))
              else sprintf("ok logLik %.6f", as.numeric(logLik(f))),
              length(w), paste(substr(unique(w), 1, 60), collapse = " ; ")))
}
cat("END\n")
