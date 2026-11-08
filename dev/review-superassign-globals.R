# Review of lane superassign: for every function in the changed files,
# compare its free variables (codetools::findGlobals) between the base
# source (main checkout, 3cd2a0c8) and the lane source. A variable that
# is free in the lane copy but not in the base copy is a reference the
# rewrite left dangling (an old local read after its definition moved
# into an environment). Also lists every local that shadows the new
# environment's name inside a nested function.
base_root <- "C:/Users/adf44/source/r/frmtmb"
lane_root <- "C:/Users/adf44/source/r/frmtmb-wt-superassign"
files <- c(
  file.path("R", c("autoscale.R", "brms-shapes.R", "ce-levels.R", "compat.R",
    "confint.R", "curve-cov.R", "fit.R", "frame.R", "generic-owners.R",
    "gp-by.R", "hazard-containers.R", "importance.R", "me.R",
    "methods-fit.R", "par-template.R", "parse.R", "predict.R", "priors.R",
    "ps.R", "re-formula.R", "simulate-new.R")),
  "extensions/frmtmb.ode/R/ode.R", "extensions/frmtmb.sample/R/loo.R",
  "extensions/frmtmb.sample/R/methods-draws.R",
  "extensions/frmtmb.sample/R/prior-draws.R",
  "extensions/frmtmb.sample/R/sample.R",
  "extensions/frmtmb.survival/R/cox-family.R")

# collect every top-level function definition of a file without
# evaluating anything but the function expressions themselves
defs <- function(path) {
  ex <- parse(path, keep.source = FALSE)
  out <- list()
  for (e in ex) {
    if (is.call(e) && as.character(e[[1L]]) %in% c("<-", "=") &&
          is.name(e[[2L]])) {
      rhs <- e[[3L]]
      nm <- as.character(e[[2L]])
      if (is.call(rhs) && identical(rhs[[1L]], as.name("function"))) {
        out[[nm]] <- eval(rhs, baseenv())
      } else if (is.call(rhs) && identical(rhs[[1L]], as.name("local"))) {
        f <- tryCatch(eval(rhs, new.env(parent = baseenv())),
                      error = function(e) NULL)
        if (is.function(f)) out[[nm]] <- f
      }
    }
  }
  out
}
glob <- function(f) {
  sort(unique(codetools::findGlobals(f, merge = FALSE)$variables))
}
cat("lib: none (source comparison)\n")
n_fun <- 0L
n_diff <- 0L
for (fl in files) {
  b <- defs(file.path(base_root, fl))
  l <- defs(file.path(lane_root, fl))
  for (nm in intersect(names(b), names(l))) {
    if (identical(deparse(b[[nm]]), deparse(l[[nm]]))) next
    n_fun <- n_fun + 1L
    gb <- glob(b[[nm]])
    gl <- glob(l[[nm]])
    new <- setdiff(gl, gb)
    gone <- setdiff(gb, gl)
    cat(sprintf("CHANGED %-28s %s\n", nm, fl))
    if (length(new)) {
      n_diff <- n_diff + 1L
      cat("  NEW FREE VARIABLES:", paste(new, collapse = ", "), "\n")
    }
    if (length(gone)) cat("  no longer free:", paste(gone, collapse = ", "), "\n")
  }
  miss <- setdiff(names(b), names(l))
  if (length(miss)) cat("  MISSING in lane:", paste(miss, collapse = ", "), "\n")
}
cat(sprintf("RESULT changed functions %d, with new free variables %d\n",
            n_fun, n_diff))
