# Lane sampleprior: which step of a seeded sampling call moves R's
# random stream, in tmbstan and rstan, and whether brms's brm(seed = )
# leaves it alone (it reaches rstan::sampling() the same way).
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
b <- deparse(body(tmbstan::tmbstan))
cat("tmbstan lines naming RNG:\n")
cat(grep("seed|runif|rnorm|sample[.(]", b, value = TRUE), sep = "\n")
m <- methods::getMethod("sampling", "stanmodel")
b2 <- deparse(m@.Data)
cat("rstan sampling lines naming RNG:\n")
cat(grep("seed|runif|rnorm|sample[.]int", b2, value = TRUE), sep = "\n")
