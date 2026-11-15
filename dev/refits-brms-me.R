# brms 2.23.0 + marginaleffects 1.0.0 on a multinomial matrix response
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), " marginaleffects",
    as.character(packageVersion("marginaleffects")), "\n")
set.seed(5)
dm <- data.frame(x = rnorm(200), N = 10L)
pr <- cbind(1, exp(0.3 + 0.5 * dm$x), exp(-0.2 * dm$x))
pr <- pr / rowSums(pr)
dm$Y <- t(sapply(1:200, function(i) rmultinom(1, 10, pr[i, ])))
colnames(dm$Y) <- c("a", "b", "c")
fb <- brm(Y | trials(N) ~ x, family = multinomial(), data = dm,
          chains = 1, iter = 400, warmup = 200, refresh = 0, seed = 1)
r <- tryCatch(marginaleffects::predictions(fb),
              warning = function(w) paste("WARN:", conditionMessage(w)),
              error = function(e) paste("ERROR:", conditionMessage(e)))
print(if (is.character(r)) r else head(as.data.frame(r)))
r2 <- tryCatch(withCallingHandlers(marginaleffects::predictions(fb),
  warning = function(w) invokeRestart("muffleWarning")),
  error = function(e) paste("ERROR:", conditionMessage(e)))
print(if (is.character(r2)) r2 else dim(r2))
str(insight::get_data(fb))
