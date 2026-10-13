# Reviewer of lane optrap: a candidate fix for the high-R^2 false
# convergence (dev/optrap-rev-r2b.R). The unit of a gaussian() or
# student() mu is taken from the sd of the residuals of the fixed design
# (mu_residuals(), which the lane already uses for skew_normal()), not
# from sd(y). Sourced after library(frmtmb); patches the namespace of
# this process only.
local({
  ns <- asNamespace("frmtmb")
  orig <- get("predictor_unit_of", ns)
  patched <- function(frame, lp, resid_cache) {
    u0 <- orig(frame, lp, resid_cache)
    if (u0 == 1 && !identical(lp[["dpar"]], "mu")) return(u0)
    fam <- frame[["spec"]]$responses[[lp[["resp"]]]]$family
    if (!identical(lp[["dpar"]], "mu") ||
          !isTRUE(fam[["family"]] %in% c("gaussian", "student"))) {
      return(u0)
    }
    # only where the lane's own rule gave a unit or would have: the same
    # exclusions (constant sigma, links, weights) are inherited from orig
    # by requiring that orig did not refuse for a structural reason
    y <- frame[["y"]][[lp[["resp"]]]]
    res <- tryCatch(mu_residuals(frame, lp[["resp"]], resid_cache),
                    error = function(e) NULL)
    if (is.null(res) || !length(res)) return(u0)
    w <- frame[["aterm_values"]][[lp[["resp"]]]][["weights"]]
    if (is.numeric(w) && length(w) == length(res) && any(w <= 0)) {
      res <- res[w > 0]
    }
    raw <- stats::sd(res)
    if (!is.finite(raw) || raw <= 0) return(u0)
    # a structural refusal of orig (se(), equate, nonlinear, link) keeps 1
    sig <- Filter(function(l) {
      identical(l[["resp"]], lp[["resp"]]) && identical(l[["dpar"]], "sigma")
    }, frame[["linpreds"]])
    if (length(sig) != 1L || !is.null(sig[[1L]][["constant"]]) ||
          !is.null(sig[[1L]][["equate"]]) ||
          !identical(lp[["link"]]$name %||% "", "identity") ||
          !is.null(lp[["constant"]]) || !is.null(lp[["equate"]]) ||
          !is.null(lp[["nl_body"]])) return(u0)
    2^round(log2(raw))
  }
  pe <- new.env(parent = ns)
  pe$orig <- orig
  environment(patched) <- pe
  utils::assignInNamespace("predictor_unit_of", patched, "frmtmb")
  cat("predictor_unit_of patched: residual sd for gaussian and student\n")
})
