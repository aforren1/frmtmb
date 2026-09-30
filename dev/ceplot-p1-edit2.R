# Lane ceplot punch 1 (m8): the method = "predict" refusal moves into
# ce_display_kind(), ahead of brms's ordinal warning.
f <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/R/conditional-effects.R"
s <- paste(readLines(f), collapse = "\n")
rep1 <- function(s, old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("expected one match, got ", n, ": ", substr(old, 1, 60))
  sub(old, new, s, fixed = TRUE)
}
s <- rep1(s, "#' did nothing, so the other layout could not be asked for at all.
#'
#' @noRd
ce_display_kind <- function(rspec, dpar, categorical) {",
"#' did nothing, so the other layout could not be asked for at all.
#'
#' `method = \"predict\"` is refused here, before brms's ordinal warning,
#' so that a call that stops does not warn first.
#'
#' @noRd
ce_display_kind <- function(rspec, dpar, categorical, method = \"epred\") {")
s <- rep1(s, "    return(\"linpred\")
  }
  if (is.null(categorical) || isTRUE(categorical)) return(\"cats\")",
"    return(\"linpred\")
  }
  if (identical(method, \"predict\")) {
    frm_stop(\"method = \\\"predict\\\" has no meaning on an ordinal family: the \",
             \"category probabilities conditional_effects() draws ARE the \",
             \"predictive distribution, so there is no further observation \",
             \"noise to add. Use method = \\\"epred\\\" (the default), or ask \",
             \"for the latent predictor with dpar = \\\"mu\\\"\", call. = FALSE)
  }
  if (is.null(categorical) || isTRUE(categorical)) return(\"cats\")")
writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], f)
cat("edited\n")
