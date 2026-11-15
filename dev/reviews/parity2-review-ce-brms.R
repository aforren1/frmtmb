# Review of lane parity2, item 3922: brms 2.23.0's default effects on the
# shapes of dev/reviews/parity2-review-data.R: get_all_effects() and the
# names of conditional_effects() on a fixed_param fit at init 0.
#   Rscript dev/reviews/parity2-review-ce-brms.R
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
options(mc.cores = 1)
cat("brms", as.character(packageVersion("brms")), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-parity2/dev/reviews/parity2-review-data.R")
try_ <- function(label, expr) {
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage"))
  if (inherits(r, "error")) cat(label, ": ERROR ", conditionMessage(r), "\n",
                                sep = "")
  invisible(r)
}
only <- commandArgs(TRUE)
if (length(only)) shapes <- shapes[only]
for (nm in names(shapes)) {
  s <- shapes[[nm]]
  cat("\n## ", nm, "\n", sep = "")
  bt <- try_(paste(nm, "brmsterms"), brmsterms(s$f, family = s$fam))
  if (!inherits(bt, "error")) {
    ae <- brms:::get_all_effects(bt, rsv_vars = brms:::rsv_vars(bt))
    cat(nm, " get_all_effects: ",
        paste(vapply(ae, paste, "", collapse = ":"), collapse = " "), "\n",
        sep = "")
  }
  fit <- try_(paste(nm, "fit"), suppressMessages(suppressWarnings(
    brm(s$f, family = s$fam, data = rd, algorithm = "fixed_param",
        chains = 1, iter = 20, warmup = 0, init = 0, refresh = 0,
        seed = 1, silent = 2))))
  if (inherits(fit, "error")) next
  ce <- try_(paste(nm, "ce"), conditional_effects(fit))
  if (!inherits(ce, "error")) {
    cat(nm, " ce: ", paste(names(ce), collapse = " "), "\n", sep = "")
  }
  if (!is.null(s$dpar)) {
    ce2 <- try_(paste(nm, "ce dpar"), conditional_effects(fit, dpar = s$dpar))
    if (!inherits(ce2, "error")) {
      cat(nm, " ce dpar=", s$dpar, ": ", paste(names(ce2), collapse = " "),
          "\n", sep = "")
    }
  }
}
