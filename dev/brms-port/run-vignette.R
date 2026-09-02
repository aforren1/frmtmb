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

ROOT <- "C:/Users/adf44/source/r/frmtmb-wt-audit"
HERE <- file.path(ROOT, "dev/brms-port")
source(file.path(HERE, "port-lib.R"))
if (MODE %in% c("spell", "v035")) source(file.path(HERE, "patches.R")) else {
  AUTO_RETRY <- list(); PATCH <- list()
}
if (MODE == "v035") {
  # Stand-ins for a gaussian default (FN-1) and lf() (FN-10). Everything
  # else is left to fail exactly as it does today.
  AUTO_RETRY <- AUTO_RETRY["default-family"]
  PATCH <- PATCH["brms_multivariate.9.1"]
}
RESDIR <- file.path(HERE, switch(MODE, spell = "results-spell",
                                 v035 = "results-v035", "results"))
dir.create(RESDIR, showWarnings = FALSE)
OUT <- file.path(RESDIR, paste0(VIG, ".rds"))
LOG <- file.path(RESDIR, paste0(VIG, ".log"))

suppressMessages(pkgload::load_all(ROOT, quiet = TRUE, export_all = FALSE))

POST_FUNS <- c(
  "summary", "plot", "conditional_effects", "conditional_smooths",
  "hypothesis", "pp_check", "predict", "fitted", "ranef", "fixef", "coef",
  "loo", "LOO", "waic", "WAIC", "add_criterion", "bayes_R2", "loo_R2",
  "marginal_effects", "marginal_smooths", "posterior_predict",
  "posterior_epred", "posterior_linpred", "as_draws_array", "nchains",
  "expose_functions", "stancode", "standata", "make_stancode", "prior_summary",
  "VarCorr", "ngrps", "launch_shinystan", "mcmc_plot", "variables",
  "residuals", "logLik", "confint", "simulate", "anova", "AIC", "BIC"
)

all_calls <- function(e, acc = character()) {
  if (!is.call(e)) return(acc)
  nm <- call_name(e)
  if (!is.na(nm)) acc <- c(acc, nm)
  for (i in seq_along(e)) if (is.call(e[[i]])) acc <- all_calls(e[[i]], acc)
  acc
}

kind_of <- function(src) {
  p <- tryCatch(parse(text = src)[[1]], error = function(e) NULL)
  if (is.null(p)) return("other")
  fns <- all_calls(p)
  if (any(fns %in% c("frm", "frm_multiple"))) return("model")
  if ("update" %in% fns) return("model")
  if (any(fns %in% POST_FUNS)) return("post")
  "other"
}

## ------------------------------------------------------------ environment
# The eval chain is  chunk env -> shim -> globalenv -> search path.
# The shim supplies ONLY data-side conveniences (brms datasets, the one
# brms data simulator, a mirror for a dead URL). No modeling or
# post-processing function is shimmed: those must resolve to frmtmb or
# fail, which is what the audit measures.
shim <- new.env(parent = globalenv())

local({
  for (d in c("kidney", "inhaler", "loss", "epilepsy")) {
    suppressWarnings(utils::data(list = d, package = "brms", envir = globalenv()))
  }
})

shim$sim_multi_mem <- tryCatch(
  getFromNamespace("sim_multi_mem", "brms"),
  error = function(e) function(nschools = 10, nstudents = 1000, change = 0.1) {
    # brms_multilevel's simulator is internal and unexported; this is the
    # documented design (two schools per student, equal weights, a share
    # of students changing school) reproduced for the audit only.
    s1 <- sample(nschools, nstudents, TRUE)
    s2 <- s1
    ch <- sample(nstudents, round(change * nstudents))
    s2[ch] <- sample(nschools, length(ch), TRUE)
    eff <- stats::rnorm(nschools, 0, 3)
    y <- 20 + 0.5 * (eff[s1] + eff[s2]) + stats::rnorm(nstudents, 0, 5)
    data.frame(s1 = s1, s2 = s2, w1 = 0.5, w2 = 0.5, y = y)
  }
)

shim$data <- function(..., package = NULL, envir = globalenv()) {
  nm <- as.character(substitute(list(...)))[-1]
  nm <- gsub('^"|"$', "", nm)
  suppressWarnings(try(utils::data(list = nm, package = package,
                                   envir = globalenv()), silent = TRUE))
  miss <- nm[!vapply(nm, exists, logical(1), envir = globalenv())]
  if (length(miss)) {
    suppressWarnings(try(utils::data(list = miss, package = "brms",
                                     envir = globalenv()), silent = TRUE))
  }
  invisible(nm)
}

# brms_multilevel points at a UCLA URL that no longer serves the file; the
# same data is the author's own mirror used by brms_distreg.
shim$read.csv <- function(file, ...) {
  if (is.character(file) && grepl("stats.idre.ucla.edu", file)) {
    file <- "https://paul-buerkner.github.io/data/fish.csv"
  }
  utils::read.csv(file, ...)
}

# Attaching brms would shadow every frmtmb generic and invalidate the run.
shim$library <- function(package, ...) {
  p <- tryCatch(as.character(substitute(package)), error = function(e) "")
  if (identical(p, "brms")) {
    message("[audit] library(brms) suppressed")
    return(invisible())
  }
  eval(bquote(base::library(.(as.name(p)))), envir = globalenv())
}
shim$require <- shim$library

## ------------------------------------------------------------------- run
env <- new.env(parent = shim)
grDevices::pdf(NULL)
on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

chunks <- extract_vignette(VIG)
res <- list()
cat("", file = LOG)

fit_summary <- function(v) {
  out <- list()
  cls <- class(v)[1]
  out$class <- cls
  if (inherits(v, "frmtmb_fit")) {
    out$loglik <- tryCatch(as.numeric(stats::logLik(v)), error = function(e) NA_real_)
    out$fixef <- tryCatch({
      f <- fixef(v)
      if (is.matrix(f)) f[, 1] else f
    }, error = function(e) NULL)
    out$conv <- tryCatch(v$opt$convergence, error = function(e) NA)
    out$sigma <- tryCatch(stats::sigma(v), error = function(e) NA_real_)
  }
  out
}

run_one <- function(src) {
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
      r <- run_one(src)
      if (r$status == "ERROR" && length(AUTO_RETRY)) {
        for (pn in names(AUTO_RETRY)) {
          alt <- tryCatch(AUTO_RETRY[[pn]](src, r$msg), error = function(e) NULL)
          if (is.null(alt) || identical(alt, src)) next
          r2 <- run_one(alt)
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
