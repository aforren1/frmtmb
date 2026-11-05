# The impsamp lane's grid: importance = against an independent reference,
# on a build chosen by the command line. One R process per chunk.
#
#   Rscript dev/impsamp-grid.R <arm> <chunk> <nchunks> [lib]
#
# Writes dev/impsamp-out/grid-<arm>-<chunk>.tsv, one row per fit. The
# arm name labels the output; the library decides the code. Same cells,
# same data seeds, same draw seed in every arm.
#
# Reference: an own adaptive Gauss-Hermite marginal likelihood per group
# (dev/impsamp-agh.R) that shares no code with the package, maximized by
# optim() from the Laplace estimates. Against quadrature = TRUE where
# that applies (one scalar block), so the reference is itself checked.
args <- commandArgs(TRUE)
arm <- args[[1]]
chunk <- as.integer(args[[2]])
nchunks <- as.integer(args[[3]])
lib <- if (length(args) >= 4L) args[[4]] else ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
source("dev/impsamp-agh.R")

cells <- impsamp_cells()
# chunk by data cell, not by row, so one process holds a cell's reference
ckey <- paste(cells$family, cells$shape, cells$size, cells$seed)
ckid <- match(ckey, unique(ckey))
cells <- cells[(ckid %% nchunks) == (chunk %% nchunks), ]
out <- sprintf("dev/impsamp-out/grid-%s-%02d.tsv", arm, chunk)
if (file.exists(out)) file.remove(out)
hdr <- TRUE
ref_cache <- new.env()
for (i in seq_len(nrow(cells))) {
  cc <- cells[i, ]
  key <- paste(cc$family, cc$shape, cc$size, cc$seed)
  d <- impsamp_data(cc$family, cc$shape, cc$size, cc$seed)
  form <- impsamp_formula(cc$family, cc$shape)
  fam <- impsamp_family(cc$family)
  if (is.null(ref_cache[[key]])) {
    lap <- tryCatch(suppressWarnings(suppressMessages(
      frm(form, family = fam, data = d))), error = function(e) e)
    qd <- if (cc$shape == "int") {
      tryCatch(suppressWarnings(suppressMessages(
        frm(form, family = fam, data = d, quadrature = TRUE))),
        error = function(e) e)
    }
    rf <- file.path("dev/impsamp-out", paste0("ref-", gsub(" ", "_", key), ".rds"))
    ref <- if (file.exists(rf)) readRDS(rf) else if (inherits(lap, "error")) NULL else {
      r0 <- tryCatch(impsamp_ref(cc$family, cc$shape, d, lap),
               error = function(e) NULL)
      if (!is.null(r0)) saveRDS(r0, rf)
      r0
    }
    ref_cache[[key]] <- list(
      lap = if (inherits(lap, "error")) NA_real_ else as.numeric(logLik(lap)),
      qd = if (is.null(qd) || inherits(qd, "error")) NA_real_ else
        as.numeric(logLik(qd)),
      ref = if (is.null(ref)) NA_real_ else ref$ll,
      refpar = ref$par)
  }
  rc <- ref_cache[[key]]
  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(
    suppressMessages(frm(form, family = fam, data = d,
                         importance = cc$N)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  res <- data.frame(family = cc$family, shape = cc$shape, size = cc$size,
                    seed = cc$seed, N = cc$N, lap = rc$lap, quad = rc$qd,
                    ref = rc$ref, secs = el)
  if (inherits(fit, "error")) {
    msg <- gsub("[\t\r\n]+", " ", conditionMessage(fit))
    res <- cbind(res, status = if (grepl("ROSE", msg)) "rose" else "error",
                 ll = NA_real_, ll_exact_at = NA_real_, mcse = NA_real_,
                 ess_min = NA_real_, rounds = NA_integer_, capped = NA,
                 moved = NA_real_, grad = NA_real_, se_na = NA_integer_,
                 d_beta0 = NA_real_, d_sd = NA_real_, rise = NA_real_, nwarn = length(warns),
                 msg = substr(msg, 1, 160))
  } else {
    im <- fit$importance
    se <- tryCatch(suppressWarnings(sqrt(diag(vcov(fit, full = TRUE)))),
                   error = function(e) NA)
    if (anyNA(se) && length(se) == 1L) {
      se <- tryCatch(suppressWarnings(summary(fit)$coefficients[, 2]),
                     error = function(e) NA)
    }
    # the same chart as the reference, read off VarCorr() and the estimates
    st <- tryCatch(impsamp_start(cc$family, cc$shape, fit), error = function(e) NULL)
    nth <- switch(cc$shape, int = 1L, slope = 3L, twoblock = 2L)
    dbeta <- if (!is.null(st) && !is.null(rc$refpar)) max(abs(st[1:2] - rc$refpar[1:2])) else NA_real_
    nsd <- min(nth, 2L)
    dsd <- if (!is.null(st) && !is.null(rc$refpar)) max(abs(st[2L + seq_len(nsd)] - rc$refpar[2L + seq_len(nsd)])) else NA_real_
    ex <- tryCatch(impsamp_exact_at(cc$family, cc$shape, d, fit),
                   error = function(e) NA_real_)
    res <- cbind(res, status = "ok", ll = as.numeric(logLik(fit)),
                 ll_exact_at = ex, mcse = im$mcse, ess_min = im$ess_min,
                 rounds = im$rounds, capped = im$capped, moved = im$moved,
                 grad = im$grad, se_na = sum(!is.finite(se)),
                 d_beta0 = dbeta, d_sd = dsd, rise = fit$opt$objective - im$start_value, nwarn = length(warns),
                 msg = substr(gsub("[\t\r\n]+", " ",
                                   paste(warns, collapse = " | ")), 1, 160))
  }
  utils::write.table(res, out, sep = "\t", row.names = FALSE, quote = TRUE,
                     col.names = hdr, append = !hdr)
  hdr <- FALSE
}
cat("DONE", arm, chunk, nrow(cells), "\n")
