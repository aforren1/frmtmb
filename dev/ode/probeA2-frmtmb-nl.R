# Probe A2: the target frmtmb spelling. Population PK where the nl body
# calls RTMBode::ode() once per subject.
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib",
            .libPaths()))
library(frmtmb)
library(RTMBode)
source("dev/ode/pk-common.R")
cat("frmtmb", as.character(packageVersion("frmtmb")), "\n")

d <- sim_pk()

form <- bf(conc ~ pk_ode(exp(lka), exp(lke), exp(lV), time, id, dose),
           lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)

cat("\n--- parse/frame dry run ---\n")
dr <- try(frm(form + gaussian(), data = d, dry_run = "frame",
              start = list(beta = c(0, log(0.25), log(8)))), silent = TRUE)
if (inherits(dr, "try-error")) {
  cat("DRY RUN FAILED:\n"); cat(as.character(dr), "\n")
} else {
  lp <- dr$linpreds[[which(vapply(dr$linpreds,
                                  function(l) !is.null(l$nl_body), TRUE))]]
  cat("nl datavars in data_list:", paste(names(lp$data_list), collapse = ", "),
      "\n")
  cat("classes:",
      paste(vapply(lp$data_list, function(x) class(x)[1], ""), collapse = ", "),
      "\n")
  cat("nl_env is globalenv:", identical(lp$nl_env, globalenv()), "\n")
}

cat("\n--- fit ---\n")
tt <- system.time(
  fit <- try(frm(form + gaussian(), data = d,
                 start = list(beta = c(0, log(0.25), log(8))),
                 se = TRUE, verbose = TRUE), silent = TRUE))
if (inherits(fit, "try-error")) {
  cat("FIT FAILED:\n"); cat(as.character(fit), "\n")
  quit(status = 0)
}
cat("elapsed:", tt[["elapsed"]], "s\n")
print(summary(fit))
saveRDS(list(fit = fit, d = d), "dev/ode/probeA2-fit.rds")
cat("PROBEA2 OK\n")
