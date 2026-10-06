# Run one transformed brms vignette against frmtmb.
#
#   Rscript run-vignette.R <vignette> [time_cap_seconds]
#
# One process per vignette on purpose: a hard crash inside a fit then
# costs one vignette's results, not the whole audit. Every expression is
# checkpointed to results/<vignette>.rds as soon as it finishes.

args <- commandArgs(trailingOnly = TRUE)
VIG <- args[1]
CAP <- if (length(args) > 1) as.numeric(args[2]) else 120
# "raw":   brm -> frm and MCMC-argument removal only.
# "spell": the same, plus the documented spelling changes in patches.R.
# "v035":  projection - only the subset of patches standing in for the
#          gaps slated to be fixed next.
MODE <- if (length(args) > 2) args[3] else "raw"

HERE <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
source(file.path(HERE, "env.R"))
source(file.path(HERE, "port-lib.R"))
if (MODE %in% c("spell", "v035", "need")) {
  source(file.path(HERE, "patches.R"))
} else {
  AUTO_RETRY <- list(); PATCH <- list()
}
# "need": only the explicit patches named in PORT_PATCHES (comma
# separated) and no retry rule. A patch on an expression whose raw
# failure was a cascade cannot be judged from the raw pass, because its
# input never existed; this pass keeps the upstream patch and drops the
# downstream one, so each is measured on its own.
if (MODE == "need") {
  keep <- strsplit(Sys.getenv("PORT_PATCHES"), ",", fixed = TRUE)[[1]]
  AUTO_RETRY <- list()
  PATCH <- PATCH[intersect(names(PATCH), keep)]
}
if (MODE == "v035") {
  # Stand-ins for a gaussian default (FN-1) and lf() (FN-10). Everything
  # else is left to fail exactly as it does today.
  AUTO_RETRY <- AUTO_RETRY["default-family"]
  PATCH <- PATCH["brms_multivariate.9.1"]
}
# "keepprior": the raw transform, except that prior = stays in the call
# and standalone prior()/set_prior() code runs. A second definition of
# the transform, measured beside the first and never in place of it:
# since 0.43.0 frm() translates a brmsprior, so a porter may keep the
# priors, and in the nonlinear vignettes the priors are what locate
# the fit.
if (MODE == "keepprior") {
  DROP_ARGS <- setdiff(DROP_ARGS, "prior")
  PRIOR_RE <- "a^"
}
RESDIR <- file.path(PORT_OUT, switch(MODE, spell = "results-spell",
                                     v035 = "results-v035",
                                     need = "results-need",
                                     keepprior = "results-keepprior",
                                     "results"))
dir.create(RESDIR, showWarnings = FALSE)
OUT <- file.path(RESDIR, paste0(VIG, ".rds"))
LOG <- file.path(RESDIR, paste0(VIG, ".log"))

# The installed build, not a source tree: the audit measures what a user
# gets from library(frmtmb). Up to 0.34.0 this was pkgload::load_all()
# on a worktree, which attaches the same exports.
suppressMessages(library(frmtmb))

## ------------------------------------------------------------ environment
# kind_of() and the shim live in port-lib.R and shim.R, which brms-fit.R
# reads too, so both runners classify and build data the same way.
source(file.path(HERE, "shim.R"))
shim <- make_shim(suppress_brms = TRUE)

## ------------------------------------------------------------------- run
env <- new.env(parent = shim)
grDevices::pdf(NULL)
on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

chunks <- extract_vignette(VIG)
res <- list()
cat("# ", port_build(), ", mode ", MODE, ", cap ", CAP, " s\n",
    sep = "", file = LOG)

fit_summary <- function(v) {
  out <- list()
  cls <- class(v)[1]
  out$class <- cls
  if (inherits(v, c("frmtmb_fit", "frmtmb_multiple"))) {
    out$loglik <- tryCatch(as.numeric(stats::logLik(v)), error = function(e) NA_real_)
    # setNames: a one-row matrix drops its row name on `[, j]`
    out$fixef <- tryCatch({
      f <- fixef(v)
      if (is.matrix(f)) stats::setNames(f[, 1], rownames(f)) else f
    }, error = function(e) NULL)
    # The estimate-plausibility table compares these to a brms fit, so
    # keep the standard errors and the variance components too, in
    # brms's own naming (fixef() and VarCorr() follow it).
    out$fixef_se <- tryCatch({
      f <- fixef(v)
      if (is.matrix(f)) stats::setNames(f[, 2], rownames(f)) else NULL
    }, error = function(e) NULL)
    # A frm_multiple() result has no fixef() method at 0.67.0; its pooled
    # table names the intercept "(Intercept)" and carries sigma's
    # log-scale row, so both are mapped to what brms's fixef() reports.
    if (inherits(v, "frmtmb_multiple") && is.null(out$fixef) &&
        !is.null(v$pooled)) {
      p <- v$pooled[!grepl("^sigma_", rownames(v$pooled)), ]
      nm <- sub("^[(]Intercept[)]$", "Intercept", rownames(p))
      out$fixef <- stats::setNames(p$estimate, nm)
      out$fixef_se <- stats::setNames(p$se, nm)
    }
    out$sds <- tryCatch({
      vc <- VarCorr(v)
      unlist(lapply(names(vc), function(g) {
        s <- vc[[g]]$sd
        if (is.null(s)) return(NULL)
        stats::setNames(s[, 1], paste0("sd_", g, "__", rownames(s)))
      }))
    }, error = function(e) NULL)
    out$conv <- tryCatch(v$opt$convergence, error = function(e) NA)
    out$sigma <- tryCatch(stats::sigma(v), error = function(e) NA_real_)
  }
  out
}

run_one <- function(src, id) {
  # one seed per expression (shim.R), so frmtmb and brms see the same data
  set.seed(port_seed(id))
  warns <- character()
  t0 <- proc.time()[["elapsed"]]
  val <- tryCatch(
    withCallingHandlers({
      setTimeLimit(elapsed = CAP, transient = TRUE)
      wv <- withVisible(eval(parse(text = src), envir = env))
      # Auto-printing matters: a print/summary method that errors is
      # exactly what a vignette reader would hit.
      if (wv$visible) utils::capture.output(print(wv$value))
      setTimeLimit()
      list(ok = TRUE, value = wv$value)
    }, warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) list(ok = FALSE, msg = conditionMessage(e))
  )
  setTimeLimit()
  secs <- round(proc.time()[["elapsed"]] - t0, 1)
  if (isTRUE(val$ok)) {
    list(status = "OK", msg = "", secs = secs, warnings = unique(warns),
         fit = fit_summary(val$value))
  } else {
    st <- if (grepl("reached elapsed time limit|reached CPU time limit", val$msg))
      "TIMEOUT" else "ERROR"
    list(status = st, msg = val$msg, secs = secs, warnings = unique(warns))
  }
}

for (k in chunks) {
  tr <- transform_code(k$code)
  for (j in seq_along(tr)) {
    t <- tr[[j]]
    id <- sprintf("%s.%d.%d", VIG, k$idx, j)
    rec <- list(id = id, vignette = VIG, chunk = k$idx, expr = j,
                header = k$header, src = t$src, xstatus = t$status,
                dropped = t$dropped, kind = kind_of(t$src))
    if (t$status %in% c("BRMS-ONLY", "PARSE-ERROR", "SETUP-SKIP")) {
      rec$status <- t$status
      rec$msg <- t$msg
      rec$secs <- 0
    } else {
      src <- t$src
      # A per-id patch is a deliberate rewrite: use it instead of the
      # mechanical transform, do not wait for a failure.
      if (!is.null(PATCH[[id]])) {
        src <- PATCH[[id]]
        rec$patch <- "explicit"
      }
      r <- run_one(src, id)
      if (r$status == "ERROR" && length(AUTO_RETRY)) {
        for (pn in names(AUTO_RETRY)) {
          alt <- tryCatch(AUTO_RETRY[[pn]](src, r$msg), error = function(e) NULL)
          if (is.null(alt) || identical(alt, src)) next
          r2 <- run_one(alt, id)
          rec$patch <- c(rec$patch, pn)
          src <- alt
          r <- r2
          if (r$status != "ERROR") break
        }
      }
      rec$src_run <- src
      rec$status <- r$status
      rec$msg <- r$msg
      rec$secs <- r$secs
      rec$warnings <- r$warnings
      rec$fit <- r$fit
    }
    res[[id]] <- rec
    cat(sprintf("[%s] %-7s %-6s %5.1fs %s %s\n", id, rec$status, rec$kind,
                rec$secs,
                if (length(rec$patch)) paste0("<", paste(rec$patch, collapse = "+"), ">") else "",
                gsub("\n", " | ", substr(rec$src_run %||% rec$src, 1, 90))),
        file = LOG, append = TRUE)
    if (nzchar(rec$msg %||% "")) {
      cat("        ! ", gsub("\n", " | ", substr(rec$msg, 1, 300)), "\n",
          sep = "", file = LOG, append = TRUE)
    }
    saveRDS(res, OUT)
  }
}
saveRDS(res, OUT)
cat("done:", VIG, length(res), "expressions\n")
