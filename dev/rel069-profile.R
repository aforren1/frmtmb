# The warning-precedence instrument of the 0.69.0 consolidation, used as
# R_PROFILE_USER by dev/rel069-prof-par.sh. One fit can be told about a
# parameter by five reports: the standard-error warning and the boundary
# message (se_report()), the separation warning (separation_check() in
# check_convergence()), the non-finite-gradient warning
# (check_convergence()) and the nonlinear flat-direction warning
# (nl_flat_message()). The rule is one report per parameter. This logs,
# per fit, the parameters each report names, in the outer parameter
# names se_report() uses, so dev/rel069-prof-sum.R can look for a
# parameter named twice.
#
# A fit is one call of fit_assembled(); its id is kept in its cache, so
# a report that waited for the first standard-error use
# (se_deferred_report()) and the check frm() runs after an autoscale
# pre-fit (se_check()) are counted with it. Lines go to
# $REL069_PROF_LOG/<pid>.txt. Logging failures are swallowed so the fit
# is unchanged.
local({
  logdir <- Sys.getenv("REL069_PROF_LOG")
  if (!nzchar(logdir)) return(invisible())
  setHook(packageEvent("frmtmb", "onLoad"), function(...) {
    ns <- asNamespace("frmtmb")
    out_file <- file.path(logdir, paste0(Sys.getpid(), ".txt"))
    st <- new.env()
    st$n <- 0L
    st$stack <- integer(0)
    tag <- Sys.getenv("REL069_PROF_FILE")
    log1 <- function(kind, pars) {
      tryCatch({
        id <- if (length(st$stack)) st$stack[length(st$stack)] else 0L
        cat(sprintf("REP\t%s\t%d\t%d\t%s\t%s\n", tag, Sys.getpid(), id,
                    kind, paste(pars, collapse = ",")),
            file = out_file, append = TRUE)
      }, error = function(e) NULL)
    }
    with_id <- function(id, expr) {
      st$stack <- c(st$stack, id)
      on.exit(st$stack <- st$stack[-length(st$stack)])
      expr
    }
    # a fit still inside fit_assembled() has no id in its cache yet; it
    # is the fit on top of the stack
    fit_id <- function(fit) {
      id <- tryCatch(fit$cache$rel069_id, error = function(e) NULL)
      if (!is.null(id)) return(id)
      if (length(st$stack)) st$stack[length(st$stack)] else 0L
    }
    o_fa <- get("fit_assembled", envir = ns)
    w_fa <- function(...) {
      st$n <- st$n + 1L
      id <- st$n
      fit <- with_id(id, {
        log1("FIT", "")
        o_fa(...)
      })
      tryCatch(if (is.environment(fit$cache)) fit$cache$rel069_id <- id,
               error = function(e) NULL)
      fit
    }
    o_sc <- get("se_check", envir = ns)
    w_sc <- function(fit, control) with_id(fit_id(fit), o_sc(fit, control))
    o_dr <- get("se_deferred_report", envir = ns)
    w_dr <- function(fit, sdr) with_id(fit_id(fit), o_dr(fit, sdr))
    # se_report()'s own split, before it reports: the SE warning takes
    # the non-boundary parameters, the boundary message the boundary
    # ones over grouping levels
    o_sr <- get("se_report", envir = ns)
    w_sr <- function(fit, lost, act) {
      tryCatch({
        ex <- get("se_explained_pars", envir = ns)(fit)
        l2 <- lost[!names(lost) %in% ex$pars &
                     !(ex$family & lost != "bound")]
        bnd <- l2 == "boundary"
        said <- bnd & get("se_par_grouped", envir = ns)(fit, names(l2))
        if (length(l2) && any(said)) log1("BOUNDARY", names(l2)[said])
        if (length(l2) && any(!bnd)) {
          log1("SE", paste0(names(l2)[!bnd], ":", l2[!bnd]))
        }
        if (length(l2) && any(bnd & !said)) {
          log1("SILENTBOUNDARY", names(l2)[bnd & !said])
        }
      }, error = function(e) log1("SEERR", conditionMessage(e)))
      o_sr(fit, lost, act)
    }
    o_sep <- get("separation_check", envir = ns)
    w_sep <- function(fit, ...) {
      s <- o_sep(fit, ...)
      if (!is.null(s)) log1("SEPARATION", s$pars)
      s
    }
    o_cc <- get("check_convergence", envir = ns)
    w_cc <- function(fit, control) {
      g <- tryCatch(drop(fit$obj$gr(fit$opt$par)), error = function(e) NULL)
      out <- o_cc(fit, control)
      tryCatch({
        if (any(grepl("gradient at the reported optimum is not finite",
                      out$warnings, fixed = TRUE)) && length(g)) {
          log1("NONFINITE", unique(names(fit$opt$par)[!is.finite(g)]))
        }
        if (length(out$warnings)) log1("CONVWARN", "")
      }, error = function(e) NULL)
      out
    }
    o_fl <- get("nl_flat_message", envir = ns)
    w_fl <- function(obj, opt, frame, fit = NULL) {
      m <- o_fl(obj, opt, frame, fit)
      if (!is.null(m)) log1("FLAT", "")
      m
    }
    reps <- list(fit_assembled = w_fa, se_check = w_sc,
                 se_deferred_report = w_dr, se_report = w_sr,
                 separation_check = w_sep, check_convergence = w_cc,
                 nl_flat_message = w_fl)
    for (nm in names(reps)) {
      f <- reps[[nm]]
      environment(f) <- environment()
      unlockBinding(nm, ns)
      assign(nm, f, envir = ns)
      lockBinding(nm, ns)
    }
  })
})
