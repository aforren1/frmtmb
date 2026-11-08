# stan_cache_write() in tests/testthat/helper-brms.R, by construction:
#   write   compile one small program through brms_stan_model() into an
#           EMPTY cache dir: the .rds appears, no .part is left
#   crash   saveRDS() that dies mid-write: no file under the final name,
#           no .part left (the old helper left a truncated .rds)
#   exists  the target already exists: replaced, no warning
#   read    (second process) read the .rds back and sample from it, the
#           cross-session load the cache exists for
#
#   Rscript dev/ci-atomicprobe.R <write|read> <cache dir> <lib1;lib2;...>
args <- commandArgs(trailingOnly = TRUE)
mode <- args[[1L]]
dir <- args[[2L]]
.libPaths(c(strsplit(args[[3L]], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = dir)
suppressPackageStartupMessages(library(testthat))
env <- new.env()
sys.source("tests/testthat/helper-brms.R", envir = env)
code <- paste(
  "data { int<lower=0> N; vector[N] y; }",
  "parameters { real mu; }",
  "model { y ~ normal(mu, 1); }", sep = "\n")
if (mode == "write") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  stopifnot(length(list.files(dir)) == 0L)
  mod <- env$brms_stan_model(code)
  cat("write: files", paste(list.files(dir), collapse = " "), "\n")
  # crash: a saveRDS() that writes some bytes and then dies
  p2 <- file.path(dir, "crash.rds")
  crash_save <- function(object, file, ...) {
    writeBin(as.raw(1:100), file)
    stop("killed mid-write")
  }
  e2 <- new.env(parent = env)
  f <- env$stan_cache_write
  environment(f) <- e2
  e2$saveRDS <- crash_save
  r <- tryCatch(f(mod, p2), error = function(e) conditionMessage(e))
  cat("crash: error '", r, "'; final name exists: ", file.exists(p2),
      "; .part files: ", length(list.files(dir, "[.]part$")), "\n", sep = "")
  # the old helper's plain write, under the same crash
  p0 <- file.path(dir, "crash-old.rds")
  tryCatch(crash_save(mod, p0), error = function(e) NULL)
  cat("crash, old write: final name exists: ", file.exists(p0),
      "; readRDS: ", inherits(try(readRDS(p0), silent = TRUE), "try-error"),
      " (TRUE = error)\n", sep = "")
  unlink(p0)
  # exists: write over a file that is already there
  p3 <- file.path(dir, "exists.rds")
  saveRDS("old", p3)
  w <- character()
  withCallingHandlers(env$stan_cache_write(mod, p3), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  cat("exists: warnings ", length(w), "; replaced: ",
      inherits(readRDS(p3), "stanmodel"), "; .part files: ",
      length(list.files(dir, "[.]part$")), "\n", sep = "")
} else {
  key <- env$brms_stan_cache_key(code)
  mod <- readRDS(file.path(dir, paste0(key, ".rds")))
  sf <- suppressWarnings(rstan::sampling(mod, data = list(N = 3L,
    y = c(-1, 0, 1)), chains = 1, iter = 20, refresh = 0, seed = 1))
  cat("read: class", class(mod), "; draws", nrow(as.matrix(sf)), "\n")
}
