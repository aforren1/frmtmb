## Probe RTMBp for the documented parallelism mechanism.
out <- file("dev/rtmbp-probe.txt", open = "wt")
p <- function(...) cat(..., "\n", file = out, sep = "")

p("R: ", R.version.string)
p("RTMB  ", as.character(packageVersion("RTMB")))
p("RTMBp ", as.character(packageVersion("RTMBp")))
p("TMB   ", as.character(packageVersion("TMB")))

## 1. exports diff
e1 <- sort(getNamespaceExports("RTMB"))
e2 <- sort(getNamespaceExports("RTMBp"))
p("\n--- exports only in RTMBp ---"); p(paste(setdiff(e2, e1), collapse = ", "))
p("--- exports only in RTMB ---");  p(paste(setdiff(e1, e2), collapse = ", "))

## 2. internal objects mentioning parallel/openmp/thread
scan_ns <- function(pkg) {
  ns <- asNamespace(pkg)
  hits <- character()
  for (f in ls(ns, all.names = TRUE)) {
    ob <- tryCatch(get(f, envir = ns), error = function(e) NULL)
    if (!is.function(ob)) next
    b <- tryCatch(paste(deparse(ob), collapse = "\n"), error = function(e) "")
    if (grepl("openmp|nthread|n_threads|parallel|omp_", b, ignore.case = TRUE))
      hits <- c(hits, f)
  }
  hits
}
p("\n--- RTMBp internals mentioning parallel/openmp/thread ---")
p(paste(scan_ns("RTMBp"), collapse = ", "))
p("--- RTMB internals mentioning parallel/openmp/thread ---")
p(paste(scan_ns("RTMB"), collapse = ", "))

## 3. registered native routines
p("\n--- RTMBp native routines matching omp/thread/parallel ---")
dl <- getLoadedDLLs()
for (nm in c("RTMB", "RTMBp")) {
  if (!nm %in% names(dl)) next
  rr <- getDLLRegisteredRoutines(dl[[nm]])
  nms <- unlist(lapply(rr, function(z) vapply(z, function(x) x$name, "")))
  p(nm, ": n routines = ", length(nms))
  p("  matches: ", paste(grep("omp|thread|parallel|config", nms, ignore.case = TRUE,
                              value = TRUE), collapse = ", "))
}

## 4. TapeConfig / MakeADFun formals in each
p("\n--- TapeConfig formals ---")
p("RTMB : ", paste(names(formals(RTMB::TapeConfig)), collapse = ", "))
p("RTMBp: ", paste(names(formals(RTMBp::TapeConfig)), collapse = ", "))
p("--- TapeConfig() current values ---")
p("RTMB : ", paste(names(RTMB::TapeConfig()), RTMB::TapeConfig(), sep = "=", collapse = " "))
p("RTMBp: ", paste(names(RTMBp::TapeConfig()), RTMBp::TapeConfig(), sep = "=", collapse = " "))
p("--- MakeADFun formals ---")
p("RTMB : ", paste(names(formals(RTMB::MakeADFun)), collapse = ", "))
p("RTMBp: ", paste(names(formals(RTMBp::MakeADFun)), collapse = ", "))

## 5. does TMB::openmp see the RTMBp dll?
p("\n--- TMB::openmp ---")
p("formals: ", paste(names(formals(TMB::openmp)), collapse = ", "))
p(paste(deparse(TMB::openmp), collapse = "\n"))

## 6. config() in each namespace?
for (nm in c("RTMB", "RTMBp")) {
  ns <- asNamespace(nm)
  p("\n--- ", nm, " has config(): ", exists("config", envir = ns, inherits = FALSE))
}

## 7. DLL: does RTMBp.dll link libgomp?
for (nm in c("RTMB", "RTMBp")) {
  f <- system.file("libs", .Platform$r_arch, paste0(nm, ".dll"), package = nm)
  p("\n", nm, " dll: ", f, "  size=", file.size(f))
}
close(out)
cat("wrote dev/rtmbp-probe.txt\n")
