## Second probe: how does RTMBp enable threads?
out <- file("dev/rtmbp-probe2.txt", open = "wt")
p <- function(...) cat(..., "\n", file = out, sep = "")

p("--- RTMBp:::.onLoad ---")
p(paste(deparse(RTMBp:::.onLoad), collapse = "\n"))
p("\n--- RTMB:::.onLoad ---")
p(if (exists(".onLoad", envir = asNamespace("RTMB"), inherits = FALSE))
    paste(deparse(get(".onLoad", envir = asNamespace("RTMB"))), collapse = "\n")
  else "RTMB has NO .onLoad (no openmp call)")

p("\n--- routines only in RTMBp ---")
library(RTMBp)
library(RTMB)
dl <- getLoadedDLLs()
rts <- function(nm) {
  rr <- getDLLRegisteredRoutines(dl[[nm]])
  sort(unlist(lapply(rr, function(z) vapply(z, function(x) x$name, ""))))
}
a <- rts("RTMB"); b <- rts("RTMBp")
p(paste(setdiff(b, a), collapse = ", "))
p("--- routines only in RTMB ---")
p(paste(setdiff(a, b), collapse = ", "))

p("\n--- TMB::config(DLL='RTMBp') ---")
cfgp <- try(TMB::config(DLL = "RTMBp"))
if (!inherits(cfgp, "try-error")) p(paste(names(cfgp), unlist(cfgp), sep = "=", collapse = "\n"))
p("\n--- TMB::config(DLL='RTMB') ---")
cfgr <- try(TMB::config(DLL = "RTMB"))
if (!inherits(cfgr, "try-error")) p(paste(names(cfgr), unlist(cfgr), sep = "=", collapse = "\n"))

p("\n--- TMB::openmp(max = TRUE) ---")
p(paste(capture.output(print(try(TMB::openmp(max = TRUE, DLL = "RTMBp")))), collapse = "\n"))
p("\n--- TMB::openmp(4, DLL='RTMBp') ---")
p(paste(capture.output(print(try(TMB::openmp(4, DLL = "RTMBp")))), collapse = "\n"))
p("\n--- TMB::openmp(4, DLL='RTMB') ---")
p(paste(capture.output(print(try(TMB::openmp(4, DLL = "RTMB")))), collapse = "\n"))

p("\n--- vignettes / doc dirs ---")
for (nm in c("RTMB", "RTMBp")) {
  p(nm, ": ", paste(list.files(system.file(package = nm), recursive = FALSE), collapse = " "))
}
close(out)
cat("done\n")
