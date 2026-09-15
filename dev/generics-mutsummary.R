# Emits the mutant matrix for item 2.5d as a block to paste verbatim.
# Every block of tests/testthat/test-scale-contract.R must appear at
# least once in the "caught by" column, or an assertion is not
# discriminating and has not been seen failing.
root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
dir <- file.path(root, "dev", "generics-out2")
labs <- c(
  M1 = "lognormal fitted() is the MEDIAN exp(mu)",
  M2 = "predict() defaults to the RESPONSE scale",
  M3 = "sigma() returns the LINK scale, log(sigma)",
  M4 = "residuals() default is formed on the LINK scale",
  M5 = "simulate() draws on the LINK scale",
  M6 = "fitted() is the median, predict(response) the mean",
  M7 = "VarCorr() is on the RESPONSE scale"
)

rd <- function(f) {
  p <- file.path(dir, paste0(f, ".txt"))
  if (!file.exists(p)) return(NULL)
  ln <- readLines(p, warn = FALSE)
  num <- function(k) {
    v <- grep(paste0("^", k, " "), ln, value = TRUE)
    if (!length(v)) return(NA_integer_)
    as.integer(trimws(sub(paste0("^", k), "", v[1])))
  }
  list(pass = num("PASS"), fail = num("FAIL"), err = num("ERROR"),
       assert = num("ASSERT"),
       blocks = trimws(sub("^BAD", "", grep("^BAD", ln, value = TRUE))))
}

cat("```\n")
cat("== item 2.5d: the scale test seen FAILING, per mutant ==\n")
cat("dev/generics-mutants.R and dev/generics-mutants2.R build each\n")
cat("mutant from this worktree, install it into its own library, and\n")
cat("dev/generics-runtests.R runs the UNMODIFIED test file against it.\n\n")
ok <- rd("fix-scale")
base <- NULL
cat(sprintf("%-4s %-50s %5s %5s\n", "mut", "what it breaks", "pass",
            "fail"))
hit <- character()
for (m in names(labs)) {
  r <- rd(paste0("mut-", m))
  if (is.null(r)) { cat(sprintf("%-4s MISSING\n", m)); next }
  cat(sprintf("%-4s %-50s %5d %5d\n", m, labs[[m]], r$pass, r$fail))
  for (b in r$blocks) cat(sprintf("       caught by: %s\n", b))
  hit <- c(hit, r$blocks)
}
cat(sprintf("\ndistinct test blocks caught by at least one mutant: %d\n",
            length(unique(hit))))
cat("```\n")
