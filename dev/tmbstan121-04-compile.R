# Can rstan compile a FRESH Stan program in this arm? Usage:
#   Rscript tmbstan121-04-compile.R <arm> <route>
# arm A puts the 2.32.10 pin on the path, arm B does not. tmbstan is
# not involved: this is rstan against StanHeaders, the pin's second
# reason. route is "rstan" (rstan::stan_model on a hand-written program)
# or "brms" (brms::brm with a program brms generates, never cached).
ROOT <- "C:/Users/adf44/source/r/tmbstan121-lib"
PIN <- "C:/Users/adf44/source/r/pinlib"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
route <- a[2]
.libPaths(if (arm == "A") c(file.path(ROOT, "A"), PIN, USER) else
            c(file.path(ROOT, "Bsrc"), USER))
cat("ARM", arm, "ROUTE", route, "rstan", format(packageVersion("rstan")),
    "StanHeaders", format(packageVersion("StanHeaders")), "from",
    dirname(find.package("StanHeaders")), "\n")
rstan::rstan_options(auto_write = FALSE)

# A per-run constant in the program text defeats any cache keyed on it.
tag <- format(as.numeric(Sys.time()) * 1000, scientific = FALSE)
t0 <- proc.time()[["elapsed"]]
res <- tryCatch({
  if (route == "rstan") {
    code <- sprintf(
      "// run %s\nparameters { real y; }\nmodel { y ~ normal(1.5, 0.25); }\n",
      tag)
    m <- rstan::stan_model(model_code = code, model_name = "t121",
                            verbose = identical(Sys.getenv("T121_VERBOSE"),
                                                "true"))
    f <- rstan::sampling(m, chains = 1, iter = 2000, seed = 7,
                         refresh = 0)
    d <- as.matrix(f)[, "y"]
    sprintf("OK mean=%.4f sd=%.4f (truth 1.5, 0.25)", mean(d), sd(d))
  } else {
    set.seed(3)
    dat <- data.frame(x = rnorm(40))
    dat$y <- 1 + 0.5 * dat$x + rnorm(40, sd = 0.7)
    fit <- brms::brm(y ~ x, data = dat, chains = 1, iter = 1000,
                     seed = 7, refresh = 0, backend = "rstan",
                     prior = brms::prior_string(sprintf("normal(0, 1%s)",
                       substr(tag, nchar(tag) - 4, nchar(tag))),
                       class = "b"))
    sprintf("OK b_x=%.4f", brms::fixef(fit)["x", "Estimate"])
  }
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("RESULT", arm, route, sprintf("%.0fs", proc.time()[["elapsed"]] - t0),
    res, "\n")
