# The fuzz spec Beta + cens + (1 | ga/gb) under nested quadrature: the
# guard's prediction against the peak working set, at several nodes.
#   Rscript dev/nestquad-betacens.R <ko> <ki> <max_gb>
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-nestquad-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-nestquad/tests/testthat"
env <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-warnings.R", "helper-fuzz.R")) {
  sys.source(file.path(wt, h), envir = env)
}
plan <- env$fuzz_plan(seed = 20260901L, size = 300L)
sp <- plan[plan$mode == "quadrature" & plan$re == "nested" &
             plan$family == "Beta", , drop = FALSE][1, ]
cat("SPEC", env$fuzz_spec_key(sp), "seed", sp$seed, "\n")
d <- env$fuzz_data(sp)
cat("rows", nrow(d$data %||% d), "\n")
dd <- if (is.data.frame(d)) d else d$data
bform <- eval(parse(text = env$fuzz_bf_text(sp)))
print(bform)
msgs <- character(0)
peak0 <- ps::ps_memory_info()[["peak_wset"]]
t0 <- proc.time()[[3]]
f <- tryCatch(withCallingHandlers(
  frm(bform, data = dd, quadrature = TRUE, verbose = TRUE,
      control = frmtmb_control(quad_nodes = as.integer(a[1:2]),
                               quad_max_gb = as.numeric(a[3]))),
  message = function(m) {
    msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
  }, warning = function(w) invokeRestart("muffleWarning")),
  error = function(e) e)
el <- proc.time()[[3]] - t0
peak1 <- ps::ps_memory_info()[["peak_wset"]]
sz <- if (!inherits(f, "error")) frmtmb:::quad_tape_size(f$obj)
cat(sprintf("BETACENS nodes %sx%s time %.1fs deltaMB %.0f tape %s %s\n",
            a[1], a[2], el, (peak1 - peak0) / 2^20,
            paste(sz, collapse = " "),
            if (inherits(f, "error")) conditionMessage(f) else
              sprintf("logLik %.8f", as.numeric(logLik(f)))))
writeLines(grep("memory|tape|done", msgs, value = TRUE))
