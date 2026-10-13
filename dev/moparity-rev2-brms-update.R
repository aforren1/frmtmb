# Reviewer, re-check: what brms 2.23.0's update(prior = ) does with a
# prior table whose rows have an empty prior, with duplicate slots, and
# with a row for a class the updated formula removes. Reads the prior
# the updated fit carries (prior_summary()); compiles only if brms needs
# to.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(11)
dd <- data.frame(x = rnorm(60), z = rnorm(60))
dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
p0 <- prior(normal(0, 0.05), class = b, coef = x) +
  prior(normal(0, 10), class = Intercept)
f0 <- brm(y ~ x + z, data = dd, prior = p0, chains = 1, iter = 400,
          refresh = 0, seed = 1)
show <- function(lab, f) {
  ps <- prior_summary(f)
  u <- ps[nzchar(ps$prior) & (ps$source == "user"), ]
  cat("\n==", lab, "\n")
  print(as.data.frame(ps)[, c("prior", "class", "coef", "source")],
        row.names = FALSE)
  if (!is.null(f$fit) && length(f$fit@sim)) {
    cat("b_x mean", format(mean(as.matrix(f)[, "b_x"]), digits = 4), "\n")
  }
}
show("stored", f0)
try1 <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) cat("\n==", lab, ": ERROR", conditionMessage(r), "\n")
  else show(lab, r)
}
# 1. the default_prior() table, rows with empty priors, and one edited
tbl <- default_prior(y ~ x + z, data = dd)
tbl$prior[tbl$class == "sigma"] <- "exponential(2)"
try1("default_prior() table, sigma edited",
     update(f0, prior = tbl, refresh = 0, seed = 1))
# 2. only the empty b/x row of that table (prior "")
bx <- tbl[tbl$class == "b" & tbl$coef == "x", ]
bx$prior <- ""
try1("one row, b x, prior empty", update(f0, prior = bx, refresh = 0,
                                        seed = 1))
# 3. duplicate slots in the new prior
try1("duplicate slots in new", update(f0, prior = prior(normal(0, 1),
  class = b, coef = x) + prior(normal(0, 2), class = b, coef = x),
  refresh = 0, seed = 1))
# 4. a row for the class the formula removes
try1("formula drops x, new row on b x",
     update(f0, formula. = ~ . - x, prior = prior(normal(0, 1), class = b,
                                                 coef = x),
            refresh = 0, seed = 1))
