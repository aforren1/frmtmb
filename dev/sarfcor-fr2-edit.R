# Fix round 2: newdata over the rows of the term's matrix (1), both
# columns of the map check (2). Part A: Rscript <this> <root>
root <- commandArgs(TRUE)[1]
edit <- function(p, pairs) {
  p <- file.path(root, p)
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}
edit("R/autocor-seam.R", list(
  # (2) both columns of the two-column call
  c("  if (is.null(Y) || !identical(dim(as.matrix(Y)), c(n, 2L)) ||
        any(abs(as.matrix(Y)[, 1L] - as.numeric(y)) >
              1e-8 * max(1, abs(as.numeric(y))))) {",
    "  # each column must be the vector call on that column alone
  y2 <- tryCatch(as.numeric(e[[\"unwhiten\"]](rev(u), numeric(n), sg, st, aux)),
                 error = function(err) NULL)
  Ym <- if (is.null(Y)) NULL else as.matrix(Y)
  if (is.null(Ym) || is.null(y2) || !identical(dim(Ym), c(n, 2L)) ||
        any(abs(Ym - cbind(as.numeric(y), y2)) >
              1e-8 * max(1, abs(as.numeric(y)), abs(y2)))) {"),
  c("  ac[[\"aux\"]] <- aux
  ac[[\"d\"]] <- n
  ac[[\"time_levels\"]] <- as.character(seq_len(n))",
    "  ac[[\"aux\"]] <- aux
  ac[[\"d\"]] <- n
  # the rows of the data the term's matrix covers, the dropped ones
  # among them: brms reads newdata of this many rows against it
  ac[[\"n_all\"]] <- n_all
  ac[[\"time_levels\"]] <- as.character(seq_len(n))"),
  c("autocor_seam_B <- function(ac, theta, mu, sigma) {
  e <- ac_block_entry(ac)",
    "autocor_seam_B <- function(ac, theta, mu, sigma) {
  theta <- autocor_seam_theta_on(ac, theta)
  e <- ac_block_entry(ac)"),
  c("autocor_seam_mean <- function(ac, theta, mu, sigma) {
  e <- ac_block_entry(ac)",
    "autocor_seam_mean <- function(ac, theta, mu, sigma) {
  theta <- autocor_seam_theta_on(ac, theta)
  e <- ac_block_entry(ac)"),
  c("autocor_seam_draw <- function(ac, theta, mu, sigma, nu = NULL) {
  e <- ac_block_entry(ac)",
    "autocor_seam_draw <- function(ac, theta, mu, sigma, nu = NULL) {
  theta <- autocor_seam_theta_on(ac, theta)
  e <- ac_block_entry(ac)"),
  c("#' Whether `newdata` can stand for a registered term's rows: the term's
#' fixed data are over the fitted rows, so `newdata` must have as many
#' rows, read in the fit's order, as brms reads `newdata` against a
#' matrix over the same number of rows. Anything else is refused, as
#' brms refuses it without a matrix over the new rows in its `newdata2`.
#'
#' @noRd
check_autocor_seam_newdata <- function(ac, newdata, what) {
  if (is.null(newdata) || NROW(newdata) == ac[[\"d\"]]) {
    return(invisible(TRUE))
  }
  frm_stop(what, \" with \", ac[[\"label\"]], \" needs newdata with the \",
           \"fitted rows: \", ac[[\"d\"]], \" rows, read as the fitted rows \",
           \"in their order (a counterfactual on the same units), where \",
           \"newdata has \", NROW(newdata), \". The term's fixed data are \",
           \"over the fitted rows. brms refuses the same call unless \",
           \"newdata2 holds the term's matrix over the new rows, which \",
           \"frmtmb does not take. frm_linpred() gives the linear \",
           \"predictor mu on any newdata, without the term\",
           call. = FALSE)
}",
    "#' The refusal of `newdata` whose rows are not the rows of the term's
#' matrix. brms 2.23.0 reads newdata against the matrix of `data2`, so
#' it takes exactly `nrow(M)` rows, the rows the fit dropped for missing
#' values among them, and refuses any other number unless `newdata2`
#' holds a matrix over the new rows (`dev/sarfcor-fr2-rows.R`).
#'
#' @noRd
stop_autocor_seam_newdata <- function(ac, m, what) {
  n_all <- ac[[\"n_all\"]] %||% ac[[\"d\"]]
  drop <- n_all - ac[[\"d\"]]
  frm_stop(what, \" with \", ac[[\"label\"]], \" needs newdata with one \",
           \"row per row of the term's matrix: \", n_all, \" rows, read as \",
           \"the rows of the data in their order (a counterfactual on the \",
           \"same units), where newdata has \", m, \".\",
           if (drop > 0L) {
             paste0(\" The fit dropped \", drop, \" of those rows for \",
                    \"missing values; newdata covers them too, as brms \",
                    \"reads it, and a row with a missing predictor gives NA.\")
           },
           \" brms refuses any other number of rows unless newdata2 holds \",
           \"the term's matrix over the new rows, which frmtmb does not \",
           \"take. frm_linpred() gives the linear predictor mu on any \",
           \"newdata, without the term\", call. = FALSE)
}

#' Whether a block covers `newdata`'s rows: a block's own rows, for a
#' draw on the rows it was built for.
#'
#' @noRd
check_autocor_seam_newdata <- function(ac, newdata, what) {
  if (is.null(newdata) || NROW(newdata) == ac[[\"d\"]]) {
    return(invisible(TRUE))
  }
  stop_autocor_seam_newdata(ac, NROW(newdata), what)
}

#' A registered block for `newdata`, read as brms reads it: one row per
#' row of the term's matrix, in order. When the fit dropped no row the
#' fitted block serves. When it dropped rows for missing values, the
#' block is built again over every row of the matrix, from the fit's
#' `data2`, and its parameters are carried over at their natural
#' values: brms's `lagsar` is rho itself, while frmtmb's internal
#' parameter maps onto an interval that depends on the rows
#' (`autocor_seam_theta_map()`). Any other number of rows is refused.
#'
#' @noRd
autocor_seam_for_newdata <- function(fit, ac, rspec, newdata, what) {
  if (is.null(newdata)) return(ac)
  n_all <- ac[[\"n_all\"]] %||% ac[[\"d\"]]
  m <- NROW(newdata)
  if (m != n_all) stop_autocor_seam_newdata(ac, m, what)
  if (n_all == ac[[\"d\"]]) return(ac)
  e <- ac_block_entry(ac)
  env <- rspec$formula_env
  data2 <- fit$data2 %||% list()
  data <- fit$data
  aux <- e[[\"build\"]](ac[[\"spec\"]], seq_len(n_all), n_all,
                      function(expr, what) {
                        lookup_structural(expr, data2, data, env, what)
                      })
  acn <- ac
  acn[[\"aux\"]] <- aux
  acn[[\"d\"]] <- n_all
  acn[[\"time_levels\"]] <- as.character(seq_len(n_all))
  acn[[\"theta_map\"]] <- autocor_seam_theta_map(e, ac, ac[[\"aux\"]], aux)
  acn
}

#' The map from a block's fitted internal parameters to those of the
#' same term built over other rows, holding the natural values fixed.
#' A term with `bounds` is mapped through them; a value outside the new
#' interval is refused (for sar(), `I - rho M` is then singular or has
#' a negative determinant over the new rows). A term with parameters
#' and no `bounds` is carried over only when its `natural` slot gives
#' the same values on both sets of rows.
#'
#' @noRd
autocor_seam_theta_map <- function(e, ac, aux0, aux1) {
  if (!ac[[\"npar\"]]) return(NULL)
  b <- e[[\"bounds\"]]
  if (!is.null(b)) {
    B0 <- b(aux0)
    B1 <- b(aux1)
    return(function(theta) {
      nat <- B0[, 1L] + (B0[, 2L] - B0[, 1L]) / (1 + exp(-theta))
      if (any(nat <= B1[, 1L] | nat >= B1[, 2L])) {
        frm_stop(\"A prediction with \", ac[[\"label\"]], \" on every row \",
                 \"of its matrix: the fitted value \",
                 paste(format(nat, digits = 4), collapse = \", \"),
                 \" is outside the interval the term has over those \",
                 \"rows, (\", paste(format(B1[, 1L], digits = 4),
                                    format(B1[, 2L], digits = 4),
                                    sep = \", \", collapse = \"; \"),
                 \")\", call. = FALSE)
      }
      stats::qlogis((nat - B1[, 1L]) / (B1[, 2L] - B1[, 1L]))
    })
  }
  function(theta) {
    n0 <- e[[\"natural\"]](theta, aux0)
    if (!isTRUE(all.equal(n0, e[[\"natural\"]](theta, aux1)))) {
      frm_stop(\"A prediction with \", ac[[\"label\"]], \" on every row of \",
               \"its matrix needs the term's parameters on those rows, \",
               \"and the term declares no bounds to carry them over by\",
               call. = FALSE)
    }
    theta
  }
}

#' A block's internal parameters as its own rows read them: the fitted
#' ones, or their image under the block's `theta_map` for a block built
#' over other rows.
#'
#' @noRd
autocor_seam_theta_on <- function(ac, theta) {
  f <- ac[[\"theta_map\"]]
  if (is.null(f)) theta else f(theta)
}")))
edit("R/predict.R", list(
  c("      check_autocor_seam_newdata(ac, newdata, \"A prediction\")
      dp <- autocor_mean_dpars(fit, rn, dpars_natural(
        fit, rspec, newdata, re_formula, allow_new_levels))
      return(response_mean(fam, dp, aterms_for_newdata(rspec, newdata)))",
    "      acn <- autocor_seam_for_newdata(fit, ac, rspec, newdata,
                                      \"A prediction\")
      dp <- dpars_natural(fit, rspec, newdata, re_formula,
                          allow_new_levels)
      dp[[\"mu\"]] <- autocor_seam_mean(acn, autocor_seam_theta(fit, ac),
                                      dp[[\"mu\"]], dp[[\"sigma\"]])
      return(response_mean(fam, dp, aterms_for_newdata(rspec, newdata)))")))
edit("R/predict-brms.R", list(
  c("      check_autocor_seam_newdata(ctx0[[\"autocor\"]], newdata, \"predict()\")",
    "      autocor_seam_for_newdata(object, ctx0[[\"autocor\"]], rspec, newdata,
                               \"predict()\")")))
edit("R/simulate-newdata.R", list(
  c("    check_autocor_seam_newdata(ac, newdata, \"A prediction or simulation\")
    return(ac)",
    "    return(autocor_seam_for_newdata(fit, ac, rspec, newdata,
                                    \"A prediction or simulation\"))")))
cat("edited", root, "\n")
