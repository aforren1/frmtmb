## Reviewer, claim 3: what brms does with the collision models of
## dev/brmsnames-rev-collide.R, from brms's own parse (no compile).
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(brms))
set.seed(12)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), x_z = rnorm(n),
                sigma_z = rnorm(n), sigma_Intercept = rnorm(n),
                Intercept = rnorm(n), y = rnorm(n), y_x = rnorm(n), y_a = rnorm(n),
                g = factor(rep(1:10, 30)), h2 = factor(rep(1:3, 100)))
d$gh2 <- factor(rep(1:6, 50))
ms <- list(
  K1 = bf(y ~ sigma_Intercept + (1 | g)),
  K2 = bf(y ~ sigma_z + (1 | g), sigma ~ z),
  K3 = bf(y ~ 1 + (1 | g:h2) + (1 | gh2)),
  K4 = mvbf(bf(y ~ x_z), bf(y_x ~ z), rescor = FALSE),
  K4b = mvbf(bf(y_a ~ x), bf(y_x ~ z), rescor = FALSE),
  K5 = bf(y ~ Intercept + x)
)
for (nm in names(ms)) {
  cat("==", nm, "==\n")
  p <- tryCatch(q(default_prior(ms[[nm]], data = d)), error = function(e) e)
  if (inherits(p, "error")) { cat("brms ERROR:", conditionMessage(p), "\n"); next }
  p <- as.data.frame(p)
  cat("rows (class coef group resp dpar):\n")
  print(p[nzchar(p$coef) | p$class %in% c("Intercept", "sigma"),
          c("class", "coef", "group", "resp", "dpar")], row.names = FALSE)
  sc <- tryCatch(q(stancode(ms[[nm]], data = d)), error = function(e) e)
  if (inherits(sc, "error")) cat("stancode ERROR:", conditionMessage(sc), "\n")
}
