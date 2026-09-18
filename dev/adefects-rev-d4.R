source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages(library(frmtmb))
TD <- file.path(tempdir(), "d4"); dir.create(TD, showWarnings = FALSE)

addr <- function(x) {
  out <- utils::capture.output(.Internal(inspect(x)))[1L]
  sub("^.*(@[0-9a-fx]+).*$", "\\1", out)
}
shared <- function(fit, what = "fit") {
  a <- addr(fit[["data"]]); b <- addr(fit$frame[["data_frame"]])
  cat(sprintf("  %-34s data@%s frame@%s shared=%s identical=%s\n",
              what, a, b, identical(a, b),
              identical(fit[["data"]], fit$frame[["data_frame"]])))
  identical(a, b)
}

mk <- function(n, seed = 20260917) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = factor(rep_len(1:8, n)))
  d$y <- rnorm(n, 1 + 0.5 * d$x, 1)
  d
}

cat("== small design, n = 40\n")
d <- mk(40)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
cat("  has data element:", !is.null(fit[["data"]]), " dim:",
    paste(dim(fit[["data"]]), collapse = "x"), "\n")
shared(fit, "fresh fit")

f1 <- file.path(TD, "fit.rds")
saveRDS(fit, f1)
fit2 <- readRDS(f1)
shared(fit2, "after saveRDS/readRDS")

lst <- list(a = fit, b = 1)
f2 <- file.path(TD, "lst.rds"); saveRDS(lst, f2)
lst2 <- readRDS(f2)
shared(lst2$a, "inside a list, saved/reloaded")

f3 <- file.path(TD, "fit.rda"); save(fit, file = f3)
e <- new.env(); load(f3, envir = e)
shared(e$fit, "after save()/load()")

up <- update(fit, formula. = y ~ x)
shared(up, "after update()")
rf <- tryCatch(refit(fit), error = function(e) e)
if (!inherits(rf, "condition")) shared(rf, "after refit()") else
  cat("  refit ERR:", substr(conditionMessage(rf), 1, 60), "\n")

bo <- tryCatch(frm_bootstrap(fit, nsim = 4, seed = 2),
               error = function(e) e)
cat("  frm_bootstrap:", if (inherits(bo, "condition"))
      substr(conditionMessage(bo), 1, 60) else class(bo)[1L], "\n")

dl <- list(d, d)
mf <- tryCatch(frm_multiple(bf(y ~ x + (1 | g)), family = gaussian(),
                            data = dl), error = function(e) e)
if (!inherits(mf, "condition")) {
  shared(mf$fits[[1L]], "inside frm_multiple fit 1")
  f4 <- file.path(TD, "mf.rds"); saveRDS(mf, f4)
  shared(readRDS(f4)$fits[[1L]], "frm_multiple after saveRDS")
} else cat("  frm_multiple ERR:",
           substr(conditionMessage(mf), 1, 60), "\n")

cat("\n== on-disk size, this arm\n")
for (n in c(40L, 20000L)) {
  dd <- mk(n)
  ff <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  p <- file.path(TD, paste0("size", n, ".rds"))
  saveRDS(ff, p)
  pu <- file.path(TD, paste0("size", n, "u.rds"))
  saveRDS(ff, pu, compress = FALSE)
  cat(sprintf("  n=%-6d gz=%9d bytes  raw=%10d bytes  objsize=%10d\n",
              n, file.size(p), file.size(pu),
              as.numeric(utils::object.size(ff))))
}

cat("\n== every other brmsfit element name\n")
bn <- c("formula", "family", "ranef", "criteria", "version", "algorithm",
        "backend", "basis", "stanvars", "model", "data", "data2", "prior",
        "fit", "exclude", "file", "save_pars", "stan_args", "threads")
for (k in bn) {
  v <- fit[[k]]
  cat(sprintf("  fit$%-12s -> %s\n", k,
              if (is.null(v)) "NULL" else paste(class(v), collapse = ",")))
}
