# Reviewer: brms 2.23.0 update() data handling, run with an empty fit.
# update.brmsfit() reads object$fit@sim for the sampler settings, so the
# empty fit gets a stub stanfit; recompile = TRUE with empty = TRUE then
# builds the updated brmsfit without compiling.
# Log: dev/drop1raw-rev-out/brms2.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages({library(brms); library(rstan)})
cat("brms", format(packageVersion("brms")), "\n")
set.seed(5)
n <- 40
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n),
                w = rnorm(n), g = factor(sample(letters[1:4], n, TRUE)))
d$z <- 1 + d$x + log(d$x2) + rnorm(n)
d$x3[c(3, 17)] <- NA
stub <- function(b) {
  b$fit <- methods::new("stanfit", sim = list(
    iter = 2, warmup = 1, chains = 1, thin = 1,
    samples = list(structure(list(), args = list(control = list())))))
  b
}
show <- function(lab, expr) {
  r <- tryCatch(suppressWarnings(expr), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-46s ERROR %s\n", lab,
                gsub("\n", " ", conditionMessage(r))))
  } else {
    cat(sprintf("%-46s nrow %d | cols %s | x2[1] %.4f\n", lab,
                nrow(r$data), paste(names(r$data), collapse = ","),
                if ("x2" %in% names(r$data)) r$data$x2[1] else NA))
  }
  invisible(r)
}
up <- function(...) update(..., recompile = TRUE, empty = TRUE)
b0 <- stub(show("fit z ~ x + log(x2) + x3, x3 NA in 2 rows",
                brm(z ~ x + log(x2) + x3, data = d, empty = TRUE)))
cat("d$x2[1] before the change:", d$x2[1], "\n")
d_orig <- d
d$x2 <- d$x2 * 10
d <- d[-(1:5), ]
show("update(b0, ~ . - x3), d changed (x2*10, -5 rows)",
     up(b0, ~ . - x3))
show("update(b0, family = student())", up(b0, family = student()))
show("update(b0, ~ . + w)  (w only in d)", up(b0, ~ . + w))
show("update(b0, ~ . + w, newdata = d)", up(b0, ~ . + w, newdata = d))
rm(d)
show("update(b0, ~ . - x) after rm(d)", up(b0, ~ . - x))
d <- d_orig
b1 <- stub(brm(z ~ poly(x2, 2) + x3, data = d, empty = TRUE))
b1u <- up(b1, ~ . - x3)
s1 <- standata(b1); s1u <- standata(b1u)
cat("brms poly(x2, 2) fit column == update column:",
    isTRUE(all.equal(unname(s1$X[, 2]), unname(s1u$X[, 2]))), "\n")
pf <- poly(d$x2, 2)[!is.na(d$x3), 1]
pc <- poly(d$x2[!is.na(d$x3)], 2)[, 1]
cat("brms fit's poly column is poly() on all 40 rows:",
    isTRUE(all.equal(unname(s1$X[, 2]), unname(pf))),
    "| on the 38 complete rows:",
    isTRUE(all.equal(unname(s1$X[, 2]), unname(pc))), "\n")
