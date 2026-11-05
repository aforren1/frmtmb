# Item 8 and 3: frm_sample(laplace = TRUE) with cores > 1 retapes on workers without the
# scope. Peak working set of this process's children, and wall. Args: arm n cores
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- args[1]; n <- as.integer(args[2]); cores <- as.integer(args[3])
rev_libs(arm)
suppressPackageStartupMessages(library(frmtmb.sample))
cat("lib sample:", dirname(system.file(package = "frmtmb.sample")), "\n")
set.seed(20261103)
x <- sort(runif(n, 0, 10)); d <- data.frame(x = x, y = sin(x) + 0.3 * rnorm(n))
fit <- suppressWarnings(frm(y ~ gp(x), d))
cat(sprintf("fit done, this process's working set %.0f MB\n", memory.size_mb <- as.numeric(system2("powershell", c("-NoProfile", "-Command", shQuote(sprintf("(Get-Process -Id %d).PeakWorkingSet64/1MB", Sys.getpid()))), stdout = TRUE)[1])))
mon <- file.path(tempdir(), "mon.ps1")
outf <- file.path(REV_OUT, sprintf("mem-%s-%d-%d.txt", arm, n, cores))
writeLines(c(sprintf("$root = %d", Sys.getpid()), "$max = 0", "while (Get-Process -Id $root -ErrorAction SilentlyContinue) {",
             "  $all = Get-CimInstance Win32_Process -Property ProcessId,ParentProcessId,WorkingSetSize", "  $kids = @($root); $grew = $true",
             "  while ($grew) { $grew = $false; foreach ($p in $all) { if ($kids -contains $p.ParentProcessId -and -not ($kids -contains $p.ProcessId)) { $kids += $p.ProcessId; $grew = $true } } }",
             "  $sum = ($all | Where-Object { $kids -contains $_.ProcessId } | Measure-Object WorkingSetSize -Sum).Sum / 1MB",
             "  if ($sum -gt $max) { $max = $sum; Set-Content -Path $env:MONOUT -Value $max }", "  Start-Sleep -Milliseconds 1500", "}"), mon)
Sys.setenv(MONOUT = outf)
system2("powershell", c("-NoProfile", "-File", shQuote(mon)), wait = FALSE)
t0 <- Sys.time()
dr <- tryCatch(suppressWarnings(suppressMessages(frm_sample(fit, laplace = TRUE, chains = 2, cores = cores, iter = 60, warmup = 30, seed = 3, refresh = 0))), error = function(e) e)
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (inherits(dr, "error")) cat("ERR", conditionMessage(dr), "\n")
Sys.sleep(3)
cat(sprintf("RME %s n %d cores %d wall %.0f s peak tree working set %s MB\n", arm, n, cores, el, paste(readLines(outf), collapse = "")))
cat("DONE\n")
