# Review of lane gpgap, item 1: where does the rewrite land? brms::stancode()
# only; no compile. For each spelling: how many nugget lines, what the edit
# changes, and that it changes nothing else.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-scope.R scope.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
e <- new.env()
sys.source("dev/equiv-specs.R", envir = local({
  x <- new.env(); x$S <- function(...) NULL; x
}), keep.source = FALSE) |> try(silent = TRUE)
rewrite <- function(code) {
  gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-6 * square(sdgp);", code,
       fixed = TRUE)
}
set.seed(1)
d <- data.frame(y = rnorm(30), y2 = rnorm(30), x = runif(30), z = runif(30),
                f = factor(rep(c("a", "b"), 15)), w = runif(30),
                ct = rpois(30, 3))
specs <- list(
  mu = list(bf(y ~ gp(x)), gaussian()),
  dpar_sigma = list(bf(y ~ 1, sigma ~ gp(x)), gaussian()),
  both_dpars = list(bf(y ~ gp(x), sigma ~ gp(z)), gaussian()),
  by_factor = list(bf(y ~ gp(x, by = f)), gaussian()),
  by_numeric = list(bf(y ~ gp(x, by = w)), gaussian()),
  gr_false = list(bf(y ~ gp(x, gr = FALSE)), gaussian()),
  twod_aniso = list(bf(y ~ gp(x, z, iso = FALSE)), gaussian()),
  scale_false = list(bf(y ~ gp(x, scale = FALSE)), gaussian()),
  approx = list(bf(y ~ gp(x, k = 10)), gaussian()),
  exact_plus_approx = list(bf(y ~ gp(x) + gp(z, k = 10)), gaussian()),
  two_terms = list(bf(y ~ gp(x) + gp(z)), gaussian()),
  mvbf = list(mvbf(bf(y ~ gp(x)), bf(y2 ~ gp(z)), rescor = FALSE), NULL),
  nl = list(bf(y ~ b1 + b2 * x, b1 ~ gp(z), b2 ~ 1, nl = TRUE), gaussian()),
  poisson = list(bf(ct ~ gp(x)), poisson()),
  matern32 = list(bf(y ~ gp(x, cov = "matern32")), gaussian()),
  matern52 = list(bf(y ~ gp(x, cov = "matern52")), gaussian()),
  exponential = list(bf(y ~ gp(x, cov = "exponential")), gaussian()),
  matern32_approx = list(bf(y ~ gp(x, cov = "matern32", k = 10)), gaussian())
)
for (nm in names(specs)) {
  s <- specs[[nm]]
  code <- tryCatch(
    if (is.null(s[[2]])) make_stancode(s[[1]], data = d) else
      make_stancode(s[[1]], family = s[[2]], data = d),
    error = function(e) paste("ERROR", conditionMessage(e)))
  if (startsWith(code, "ERROR")) {
    cat(sprintf("%-20s %s\n", nm, substr(code, 1, 120)))
    next
  }
  L <- strsplit(code, "\n")[[1]]
  new <- strsplit(rewrite(code), "\n")[[1]]
  ch <- which(L != new)
  fun <- vapply(ch, function(i) {
    h <- rev(grep("^\\s*(vector|matrix|real) [a-z_]+\\(", L[seq_len(i)], value = TRUE))
    if (length(h)) trimws(sub("\\(.*", "", h[1])) else "?"
  }, "")
  cat(sprintf(paste0("%-20s lines changed %d, in %s; all 12 nugget spellings ",
                     "left after edit: %d; sdgp is in scope: %s; fn calls: %s\n"),
              nm, length(ch), paste(fun, collapse = ","),
              sum(grepl("1e-12", new, fixed = TRUE)),
              if (length(ch)) all(vapply(ch, function(i) {
                any(grepl("real sdgp", L[seq_len(i)], fixed = TRUE))
              }, NA)) else NA,
              paste(unique(sub(".*(gp_[a-z_0-9]+|gpa)\\(.*", "\\1",
                               grep("= *(gp_[a-z_0-9]+|gpa)\\(|\\+= *(gp_[a-z_0-9]+|gpa)\\(", L, value = TRUE))),
                    collapse = ",")))
}
code <- make_stancode(bf(y ~ gp(x)), family = gaussian(), data = d)
L <- strsplit(code, "\n")[[1]]
i <- grep("cov[n, n] += 1e-12", L, fixed = TRUE)
cat(L[(i - 16):(i + 4)], sep = "\n")
