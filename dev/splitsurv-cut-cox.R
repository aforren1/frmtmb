# Cut cox() out of core's R/families.R, saving what was cut so that
# frmtmb.survival's R/cox.R is built from the same text. Run once.
f <- "R/families.R"
x <- readLines(f)
one <- function(pat, fixed = TRUE) {
  i <- grep(pat, x, fixed = fixed)
  if (length(i) != 1L) stop("expected one line for: ", pat)
  i
}
cut <- list()
take <- function(from, to, label) {
  cut[[label]] <<- x[from:to]
  x[from:to] <<- NA_character_
}
# the family, its bases and its cached designs
a <- one("# --- Cox proportional hazards ---")
b <- one("#' Call a family's log-CDF, passing the family-level extra parameters to")
take(a, b - 1L, "family")
# cox_baseline()
a <- one("#' The fitted baseline-hazard simplex of a `cox()` fit.")
b <- one("#' The families whose response is an `n x K` matrix of counts over")
take(a, b - 1L, "baseline")
# the registry entry
take(one("  cox                       = fam_cox,"),
     one("  cox                       = fam_cox,"), "registry")
# the doc section
a <- one("#' @section Cox proportional hazards:")
b <- one("#' @section Quantile regression inference:")
take(a, b - 1L, "docs")
# the example
a <- one("#' # proportional hazards with a spline baseline; (1 | g) is a frailty")
b <- one("#' cox_baseline(cox_fit)              # the baseline hazard weights")
# the blank example line before it goes too
take(a - 1L, b, "example")
# the constructor and its parameters
a <- one("#' @param df For `cox()`: the number of M-spline basis functions in the")
b <- one("cox <- function(link = \"log\", df = 5, degree = 3, intercept = TRUE) {")
e <- b + which(x[(b + 1L):length(x)] == "}")[1L]
take(a - 1L, e + 1L, "ctor")
x <- x[!is.na(x)]
writeLines(x, f)
dir.create("dev/splitsurv-out", showWarnings = FALSE)
con <- file("dev/splitsurv-out/cox-cut.R", "w")
for (nm in names(cut)) {
  writeLines(c(paste0("##### ", nm), cut[[nm]]), con)
}
close(con)
cat(vapply(cut, length, 0L), "\n")
