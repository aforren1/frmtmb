LIB <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
say <- function(...) cat(sprintf(...))
q <- function(...) suppressMessages(suppressWarnings(...))
q(library(frmtmb))
m <- cbind(a = rnorm(50))
say("== the break: a caller who named the first argument ==\n")
for (call in c('posterior_summary(object = m)', 'posterior_summary(x = m)',
               'nvariables(x = m)', 'VarCorr(x = 1)')) {
  r <- tryCatch({ eval(parse(text = call)); "OK" },
                error = function(e) paste("ERR:", substr(conditionMessage(e), 1, 60)))
  say("  %-32s %s\n", call, r)
}
say("\n== where S3method(VarCorr, frmtmb_fit) actually landed ==\n")
for (p in c("nlme", "generics", "frmtmb")) {
  tb <- tryCatch(get(".__S3MethodsTable__.", envir = asNamespace(p)),
                 error = function(e) NULL)
  if (is.null(tb)) next
  hits <- grep("frmtmb", ls(tb), value = TRUE)
  say("  %-10s : %s\n", p, paste(hits, collapse = ", "))
}
say("\n== adoption loop cost, independently ==\n")
fa <- get("frm_adopt_generics", envir = asNamespace("frmtmb"))
hk <- function() sum(vapply(c("brms","lme4","posterior","loo","rstantools",
  "bayesplot"), function(p) length(getHook(packageEvent(p,"onLoad"))), 0L))
say("  hooks after load: %d\n", hk())
invisible(fa("frmtmb")); h1 <- hk()
for (i in 1:3) invisible(fa("frmtmb"))
h4 <- hk()
say("  hooks after 1 more call: %d ; after 4: %d  (must be equal)\n", h1, h4)
n <- 1024L
t <- system.time(for (i in seq_len(n)) fa("frmtmb"))[["elapsed"]]
say("  %d reps %.3f s  %.1f us each\n", n, t, 1e6 * t / n)
cat("GENREVDONE\n")
