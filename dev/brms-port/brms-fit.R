# Fit one vignette's models with brms itself, for the estimate
# comparison.
#
#   Rscript brms-fit.R <vignette> [iter] [chains]
#
# Why: "Estimate plausibility" in dev/brms-vignette-port.md compared
# frmtmb's maximum-likelihood estimates with the posterior means the
# vignettes PRINT. Several vignettes simulate their data without a
# seed, so the printed numbers belong to data nobody can rebuild. This
# runner evaluates the vignette's ORIGINAL code with brms attached,
# under the per-expression seeds that run-vignette.R uses, so both
# sides fit the same rows, and records each brmsfit's summaries.
#
# The vignette code runs as written except for the sampler budget:
# chains, iter, warmup and cores are set on every brm(), brm_multiple()
# and update() call, file = is dropped, and seed = is added where the
# vignette gives none. Priors, control and everything else stay.
# Post-processing expressions are skipped: only the fits are needed.
#
# Output: <PORT_OUT>/brms/<vignette>.rds and .log.
HERE <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
source(file.path(HERE, "env.R"))
source(file.path(HERE, "port-lib.R"))
source(file.path(HERE, "shim.R"))
args <- commandArgs(trailingOnly = TRUE)
VIG <- args[1]
ITER <- if (length(args) > 1) as.integer(args[2]) else 1000L
CHAINS <- if (length(args) > 2) as.integer(args[3]) else 2L
CAP <- 3600
# iter = 0 keeps the vignette's own sampler settings, and brms's
# defaults (4 chains of 2000) where the vignette gives none: the budget
# for a model the cheap budget does not converge on. Its records go to
# brms-vig/, which plausibility.R prefers where both exist.
VIGBUDGET <- identical(ITER, 0L)
RESDIR <- file.path(PORT_OUT, if (VIGBUDGET) "brms-vig" else "brms")
dir.create(RESDIR, showWarnings = FALSE)
OUT <- file.path(RESDIR, paste0(VIG, ".rds"))
LOG <- file.path(RESDIR, paste0(VIG, ".log"))

suppressMessages(library(brms))
options(mc.cores = CHAINS, brms.backend = "rstan")
cat(if (VIGBUDGET) sprintf(
  "# brms %s, rstan %s; the vignette's own sampler settings\n",
  utils::packageVersion("brms"), utils::packageVersion("rstan")) else
    sprintf("# brms %s, rstan %s; %d chains x %d iter (%d warmup)\n",
            utils::packageVersion("brms"), utils::packageVersion("rstan"),
            CHAINS, ITER, ITER %/% 2L), file = LOG)

budget <- function(e, id) {
  if (!is.call(e)) return(e)
  nm <- call_name(e)
  if (!is.na(nm) && nm %in% c("brm", "brm_multiple", "update")) {
    # `[[<-` with NULL on an absent name of a call is an error, so drop
    # by subsetting
    if (VIGBUDGET) {
      ch <- if (is.null(e$chains)) 4L else e$chains
      if (!is.null(names(e))) {
        e <- e[!names(e) %in% c("file", "file_refit", "cores", "refresh")]
      }
      e$cores <- ch
    } else {
      if (!is.null(names(e))) {
        e <- e[!names(e) %in% c("file", "file_refit", "iter", "warmup",
                                "chains", "cores", "refresh")]
      }
      e$chains <- CHAINS
      e$iter <- ITER
      e$warmup <- ITER %/% 2L
      e$cores <- CHAINS
    }
    e$refresh <- 0
    if (is.null(e$seed)) e$seed <- port_seed(id)
  }
  for (i in seq_along(e)) {
    if (!is.null(e[[i]]) && is.call(e[[i]])) e[[i]] <- budget(e[[i]], id)
  }
  e
}

brms_summary <- function(v) {
  out <- list(class = class(v)[1])
  fx <- tryCatch(brms::fixef(v), error = function(e) NULL)
  if (!is.null(fx)) {
    # setNames: a one-row matrix drops its row name on `[, j]`
    out$fixef <- stats::setNames(fx[, "Estimate"], rownames(fx))
    out$fixef_sd <- stats::setNames(fx[, "Est.Error"], rownames(fx))
  }
  out$sds <- tryCatch({
    vc <- brms::VarCorr(v)
    unlist(lapply(names(vc), function(g) {
      s <- vc[[g]]$sd
      if (is.null(s)) return(NULL)
      stats::setNames(s[, "Estimate"], paste0("sd_", g, "__", rownames(s)))
    }))
  }, error = function(e) NULL)
  out$sds_sd <- tryCatch({
    vc <- brms::VarCorr(v)
    unlist(lapply(names(vc), function(g) {
      s <- vc[[g]]$sd
      if (is.null(s)) return(NULL)
      stats::setNames(s[, "Est.Error"], paste0("sd_", g, "__", rownames(s)))
    }))
  }, error = function(e) NULL)
  out$sigma <- tryCatch({
    ps <- brms::posterior_summary(v, variable = "^sigma", regex = TRUE)
    ps[, "Estimate"]
  }, error = function(e) NULL)
  out$rhat_max <- tryCatch(max(brms::rhat(v), na.rm = TRUE),
                           error = function(e) NA_real_)
  out$divergent <- tryCatch({
    np <- brms::nuts_params(v)
    sum(np$Value[np$Parameter == "divergent__"])
  }, error = function(e) NA_real_)
  out$ndraws <- tryCatch(brms::ndraws(v), error = function(e) NA_integer_)
  out
}

shim <- make_shim(suppress_brms = FALSE)
env <- new.env(parent = shim)
grDevices::pdf(NULL)
res <- list()
for (k in extract_vignette(VIG)) {
  tr <- transform_code(k$code)
  ex <- parse(text = paste(k$code, collapse = "\n"), keep.source = FALSE)
  for (j in seq_along(tr)) {
    id <- sprintf("%s.%d.%d", VIG, k$idx, j)
    kind <- kind_of(tr[[j]]$src)
    rec <- list(id = id, kind = kind)
    if (tr[[j]]$status %in% c("SETUP-SKIP", "PARSE-ERROR") ||
        kind == "post") {
      rec$status <- "SKIPPED"
      res[[id]] <- rec
      next
    }
    e <- budget(ex[[j]], id)
    set.seed(port_seed(id))
    t0 <- proc.time()[["elapsed"]]
    warns <- character()
    val <- tryCatch(withCallingHandlers({
      setTimeLimit(elapsed = CAP, transient = TRUE)
      v <- eval(e, envir = env)
      setTimeLimit()
      list(ok = TRUE, value = v)
    }, warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(err) list(ok = FALSE, msg = conditionMessage(err)))
    setTimeLimit()
    rec$secs <- round(proc.time()[["elapsed"]] - t0, 1)
    rec$status <- if (isTRUE(val$ok)) "OK" else "ERROR"
    rec$msg <- val$msg %||% ""
    rec$warnings <- unique(warns)
    if (isTRUE(val$ok) && inherits(val$value, "brmsfit")) {
      rec$fit <- brms_summary(val$value)
    }
    res[[id]] <- rec
    cat(sprintf("[%s] %-6s %-6s %7.1fs %s\n", id, rec$status, kind,
                rec$secs, gsub("\n", " ", substr(deparse1(ex[[j]]), 1, 90))),
        file = LOG, append = TRUE)
    if (nzchar(rec$msg)) cat("        ! ", substr(rec$msg, 1, 300), "\n",
                             sep = "", file = LOG, append = TRUE)
    if (!is.null(rec$fit)) {
      cat(sprintf("        rhat_max %.3f  divergent %s  ndraws %s\n",
                  rec$fit$rhat_max, rec$fit$divergent, rec$fit$ndraws),
          file = LOG, append = TRUE)
    }
    saveRDS(res, OUT)
  }
}
saveRDS(res, OUT)
cat("done:", VIG, length(res), "expressions\n", file = LOG, append = TRUE)
