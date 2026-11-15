# Lane parity2, fix round: brms 2.23.0 on (1) the names of a smooth of a
# non-syntactic column, s(`my x`); (2) the display order of y ~ f * x;
# (3) cumulative() with center = FALSE (the stanc failure filed as
# brms-52). fixed_param fits at init 0, 20 draws, seed 1.
#   Rscript dev/parity2-brms-probe4.R > dev/parity2-out/brms-probe4.txt 2>&1
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
options(mc.cores = 1)
cat("brms", as.character(packageVersion("brms")), "\n")
bfix <- function(formula, data, ...) {
  suppressMessages(suppressWarnings(
    brm(formula, data = data, algorithm = "fixed_param", chains = 1,
        iter = 20, warmup = 0, init = 0, refresh = 0, seed = 1,
        silent = 2, ...)))
}
set.seed(4)
d <- data.frame(x = runif(80, 0.1, 2), f = factor(rep(c("a", "b"), 40)))
d$y <- sin(2 * d$x) + (d$f == "b") + rnorm(80, sd = 0.3)
d[["my x"]] <- d$x
d$yo <- cut(d$y, quantile(d$y, c(0, 0.3, 0.6, 1)), include.lowest = TRUE,
            labels = FALSE)

cat("\n#### s(`my x`)\n")
r <- tryCatch(bfix(y ~ s(`my x`), d), error = function(e) e)
if (inherits(r, "error")) cat("ERROR", conditionMessage(r), "\n") else {
  cat("variables:", paste(variables(r), collapse = " "), "\n")
}
pr <- tryCatch(get_prior(y ~ s(`my x`), d), error = function(e) e)
if (!inherits(pr, "error")) print(pr[pr$class == "sds", 1:4])

cat("\n#### y ~ f * x default displays\n")
fx <- bfix(y ~ f * x, d)
cat("ce:", paste(names(suppressWarnings(conditional_effects(fx))),
                 collapse = " "), "\n")
cat("ce effects = 'f:x':",
    paste(names(conditional_effects(fx, effects = "f:x")), collapse = " "),
    "\n")

cat("\n#### cumulative() with center = FALSE\n")
r <- tryCatch(bfix(bf(yo ~ x, center = FALSE), d, family = cumulative()),
              error = function(e) e)
cat(if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else
  "fits", "\n")
cat(grep("ordered_logistic", strsplit(stancode(bf(yo ~ x, center = FALSE),
                                               d, family = cumulative()),
                                      "\n")[[1]], value = TRUE), sep = "\n")
cat("done\n")
