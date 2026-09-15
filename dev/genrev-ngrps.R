.libPaths(c("C:/Users/adf44/source/r/genrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
say <- function(...) cat(sprintf(...))
q <- function(...) suppressMessages(suppressWarnings(...))
q(loadNamespace("lme4")); q(loadNamespace("brms"))
tb <- function(p) get(".__S3MethodsTable__.", envir = asNamespace(p))
say("lme4::ngrps and brms::ngrps identical: %s\n",
    identical(lme4::ngrps, brms::ngrps))
say("ngrps.brmsfit in brms table: %s ; in lme4 table: %s\n",
    exists("ngrps.brmsfit", envir = tb("brms"), inherits = FALSE),
    exists("ngrps.brmsfit", envir = tb("lme4"), inherits = FALSE))
say("ngrps.merMod in lme4 table/ns: %s ; in brms table: %s\n",
    exists("ngrps.merMod", envir = tb("lme4"), inherits = FALSE) ||
      exists("ngrps.merMod", envir = asNamespace("lme4"), inherits = FALSE),
    exists("ngrps.brmsfit", envir = tb("brms"), inherits = FALSE) &&
      exists("ngrps.merMod", envir = tb("brms"), inherits = FALSE))
say("so a single import could repair: %s\n",
    "whichever owner's table it lands in, never both")
say("\n== refit: same question ==\n")
say("lme4::refit and generics::refit identical: %s\n",
    identical(lme4::refit, generics::refit))
say("refit.merMod in lme4: %s ; in generics: %s\n",
    exists("refit.merMod", envir = tb("lme4"), inherits = FALSE) ||
      exists("refit.merMod", envir = asNamespace("lme4"), inherits = FALSE),
    exists("refit.merMod", envir = tb("generics"), inherits = FALSE))
say("\n== lme4 recursive dependency count ==\n")
d <- tools::package_dependencies("lme4", db = installed.packages(),
                                 recursive = TRUE)[[1]]
say("lme4 recursive deps: %d -> %s\n", length(d), paste(sort(d), collapse=","))
n <- tools::package_dependencies("nlme", db = installed.packages(),
                                 recursive = TRUE)[[1]]
say("nlme recursive deps: %d -> %s\n", length(n), paste(sort(n), collapse=","))
g <- tools::package_dependencies("generics", db = installed.packages(),
                                 recursive = TRUE)[[1]]
say("generics recursive deps: %d -> %s\n", length(g), paste(sort(g), collapse=","))
cat("GENREVDONE\n")
