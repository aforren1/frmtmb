# Reviewer of the 0.69.0 consolidation, used as R_PROFILE_USER: number
# fits as dev/rel069-profile.R numbers them (one per fit_assembled()
# call) and, for the fit ids listed in RELREV_TARGETS (comma separated),
# check independently that each report the precedence run logged is
# right: a parameter in the boundary message loses at most 1e-6 of
# log-likelihood when moved toward its edge with every other parameter
# held (an upper bound on the profile loss), and a parameter in the
# standard-error warning loads a direction of the outer Hessian (taken
# here by optimHess() on obj$fn and obj$gr) that the likelihood does
# not curve along. Output to $RELREV_OUT.
local({
  out <- Sys.getenv("RELREV_OUT")
  tg <- as.integer(strsplit(Sys.getenv("RELREV_TARGETS"), ",")[[1]])
  if (!nzchar(out) || !length(tg)) return(invisible())
  setHook(packageEvent("frmtmb", "onLoad"), function(...) {
    ns <- asNamespace("frmtmb")
    st <- new.env(); st$n <- 0L
    say <- function(...) cat(..., "\n", sep = "", file = out, append = TRUE)
    analyse <- function(fit, id) {
      say("== fit ", id, ": ", paste(deparse(fit$call$formula %||%
                                               fit$formula), collapse = " "))
      an <- tryCatch(get("se_fit_analysis", ns)(fit, force = TRUE),
                     error = function(e) NULL)
      lost <- an$lost
      say("lost: ", paste(names(lost), lost, sep = "=", collapse = " "))
      ex <- tryCatch(get("se_explained_pars", ns)(fit),
                     error = function(e) NULL)
      say("explained by another report: ", paste(ex$pars, collapse = " "),
          " family=", isTRUE(ex$family))
      obj <- fit$obj; p0 <- fit$opt$par
      on <- get("outer_par_names", ns)(fit)
      f0 <- obj$fn(p0)
      say("code ", fit$opt$convergence, " objective ", format(f0, digits = 12),
          " max|grad| ", format(max(abs(obj$gr(p0))), digits = 3))
      for (nm in names(lost)[lost == "boundary"]) {
        j <- match(nm, on)
        r <- vapply(c(-log(7), -3, -10, log(7), 1, -1), function(s) {
          p <- p0; p[j] <- p[j] + s; obj$fn(p) - f0
        }, 0)
        say(sprintf(paste0("  boundary %s = %.4g: loss moving it by",
                           " -log7 %.2e, -3 %.2e, -10 %.2e | +log7 %.2e,",
                           " +1 %.2e, -1 %.2e"), nm, p0[j], r[1], r[2],
                    r[3], r[4], r[5], r[6]))
      }
      se_nm <- names(lost)[lost != "boundary"]
      mt <- tryCatch(get("summary_mo_frame", ns)(fit, 0.95), error = function(e) NULL)
      if (!is.null(mt)) for (r in seq_len(nrow(mt))) say("  mo weight ", rownames(mt)[r], " = ", format(mt[r, 1], digits = 4), " se ", format(mt[r, 2], digits = 3))
      # one coordinate at a time, the others held: an upper bound on the
      # profile loss, and the gradient there
      g0 <- tryCatch(obj$gr(p0), error = function(e) rep(NA, length(p0)))
      for (nm in se_nm) {
        j <- match(nm, on)
        r <- vapply(c(-1, -0.1, 0.1, 1), function(s) {
          p <- p0; p[j] <- p[j] + s
          tryCatch(obj$fn(p) - f0, error = function(e) NA_real_)
        }, 0)
        say(sprintf(paste0("  SE-lost %s (%s) = %.4g, gradient %.3g: loss",
                           " moving it by -1 %.2e, -0.1 %.2e, +0.1 %.2e,",
                           " +1 %.2e"), nm, lost[[nm]], p0[j], g0[j],
                    r[1], r[2], r[3], r[4]))
      }
      # the outer Hessian over every parameter but the boundary and
      # non-finite ones (held at the estimate): the SE report is about
      # the others, and a step of a log sd at its edge can be infinite
      held <- names(lost)[lost %in% c("boundary", "nonfinite")]
      idx <- which(!on %in% held)
      if (!all(is.finite(g0[idx]))) say("  gradient not finite at: ", paste(on[idx][!is.finite(g0[idx])], collapse = " "))
      if (length(se_nm) && all(is.finite(g0[idx]))) {
        fs <- function(q) { p <- p0; p[idx] <- q; obj$fn(p) }
        gs <- function(q) { p <- p0; p[idx] <- q; obj$gr(p)[idx] }
        H <- tryCatch(stats::optimHess(p0[idx], fs, gs),
                      error = function(e) { say("  optimHess failed: ", conditionMessage(e)); NULL })
        if (!is.null(H)) {
          D <- sqrt(abs(diag(H))); D[D == 0] <- 1
          e <- eigen((H + t(H)) / 2 / outer(D, D), symmetric = TRUE)
          k <- length(e$values)
          on2 <- on[idx]
          say("  unit-diagonal Hessian over ", k, " parameters (held: ",
              paste(held, collapse = " "), "); smallest 4 / largest: ",
              paste(format(rev(e$values)[1:min(4, k)] / max(abs(e$values)),
                           digits = 3), collapse = " "))
          for (i in rev(order(-e$values))[1:min(2, k)]) {
            v <- e$vectors[, i]
            big <- order(-abs(v))[1:min(5, k)]
            say("   direction of eigenvalue ",
                format(e$values[i] / max(abs(e$values)), digits = 3),
                " loads: ", paste0(on2[big], "=", format(v[big], digits = 3),
                                   collapse = " "))
            dv <- v / D
            dv <- dv / max(abs(dv))
            r <- vapply(c(0.01, 0.1, 1), function(s) fs(p0[idx] + s * dv) - f0, 0)
            say("   loss along it at steps 0.01, 0.1, 1 (largest move): ",
                paste(format(r, digits = 3), collapse = " "))
          }
          fl <- abs(e$values) <= 1e-6 * max(abs(e$values))
          for (nm in intersect(se_nm, on2)) {
            j <- match(nm, on2)
            say(sprintf("  SE-lost %s (%s): projection onto the subspace below 1e-6 %.3g",
                        nm, lost[[nm]], sqrt(sum(e$vectors[j, fl]^2))))
          }
        }
      }
      obj$fn(p0)
      invisible()
    }
    o_fa <- get("fit_assembled", envir = ns)
    w_fa <- function(...) {
      st$n <- st$n + 1L
      id <- st$n
      fit <- o_fa(...)
      if (id %in% tg) tryCatch(analyse(fit, id), error = function(e)
        say("analysis error at fit ", id, ": ", conditionMessage(e)))
      fit
    }
    environment(w_fa) <- environment()
    unlockBinding("fit_assembled", ns)
    assign("fit_assembled", w_fa, envir = ns)
    lockBinding("fit_assembled", ns)
  })
})
