# Lane fixes, punch round 1: the pre-punch guard, check_nl_sum_identified()
# as it was delivered (symbolic stats::D() test), kept so the new tests
# can be run against it. dev/fixes-p0-run-test.R injects it into the
# lane build's namespace in place of check_nl_identified().
p0_check_nl_sum_identified <- function(spec, frame, prior) {
  lps <- frame[["linpreds"]] %||% list()
  bodies <- Filter(function(lp) !is.null(lp[["nl_body"]]), lps)
  if (!length(bodies)) return(invisible(NULL))
  reads <- list()
  for (lp in bodies) {
    for (v in all.vars(lp[["nl_body"]])) {
      key <- linpred_key(lp[["resp"]], v)
      reads[[key]] <- (reads[[key]] %||% 0L) + 1L
    }
  }
  prior_idx <- NULL
  for (lp in bodies) {
    r <- lp[["resp"]]
    pars <- intersect(lp[["nl_pars"]] %||% spec$responses[[r]]$nlpars,
                      all.vars(lp[["nl_body"]]))
    ok <- vapply(pars, function(p) {
      lpp <- lps[[linpred_key(r, p)]]
      !is.null(lpp) && is.null(lpp[["constant"]]) &&
        identical(lpp[["par"]], "beta") && length(lpp[["idx"]]) > 0L &&
        identical(reads[[linpred_key(r, p)]], 1L)
    }, NA)
    pars <- pars[ok]
    if (length(pars) < 2L) next
    der <- lapply(pars, function(p) {
      tryCatch(stats::D(lp[["nl_body"]], p), error = function(e) NULL)
    })
    keep <- !vapply(der, is.null, NA)
    pars <- pars[keep]
    der <- der[keep]
    while (length(pars) >= 2L) {
      same <- vapply(der, identical, NA, der[[1L]])
      grp <- pars[same]
      pars <- pars[!same]
      der <- der[!same]
      if (length(grp) < 2L) next
      X <- do.call(cbind, lapply(grp, function(p) {
        as.matrix(lps[[linpred_key(r, p)]][["X"]])
      }))
      idx <- unlist(lapply(grp, function(p) lps[[linpred_key(r, p)]][["idx"]]))
      if (ncol(X) != length(idx)) next
      q <- qr(X, tol = 1e-7)
      if (q$rank >= ncol(X)) next
      if (is.null(prior_idx)) {
        prior_idx <- if (is.null(prior)) integer(0) else {
          ent <- resolve_prior_input(list(frame = frame, spec = spec),
                                     prior)$entries
          unlist(lapply(ent, function(e) {
            if (identical(e$comp, "beta")) e$idx
          }))
        }
      }
      pen <- idx %in% prior_idx
      if (any(pen)) {
        sv <- svd(X)
        null <- sv$v[, seq.int(q$rank + 1L, ncol(X)), drop = FALSE]
        if (qr(null[pen, , drop = FALSE], tol = 1e-7)$rank == ncol(null)) {
          next
        }
      }
      labels <- unlist(lapply(grp, function(p) {
        cn <- colnames(lps[[linpred_key(r, p)]][["X"]])
        paste0(p, "_", sub("(Intercept)", "Intercept", cn, fixed = TRUE))
      }))
      alias <- labels[q$pivot[seq.int(q$rank + 1L, ncol(X))]]
      frm_stop("The nonlinear parameters ", paste0("'", grp, "'",
                                                  collapse = " and "),
               " enter `", deparse1(lp[["nl_body"]]), "` only through ",
               "their sum, and their formulas share terms, so ",
               paste(alias, collapse = ", "), " cannot be told apart ",
               "from the other parameters' coefficients: every split of ",
               "the sum has the same likelihood, and no starting value ",
               "fits this model. brms samples it only when priors ",
               "identify it. Put a prior on the coefficients of one of ",
               "them, e.g. prior = set_prior(\"normal(0, 1)\", nlpar = \"",
               grp[length(grp)], "\"), or drop the shared terms from all ",
               "but one of their formulas", call. = FALSE)
    }
  }
  invisible(NULL)
}
