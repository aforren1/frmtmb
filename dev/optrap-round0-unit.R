# Lane optrap, punch round 1: the round-0 unit (the family's own start
# of sigma, sd(y) for gaussian()), patched into this process's frmtmb
# namespace, so a test of the high-R^2 design can be seen to fail on the
# round-0 rule with everything else of the final build. Sourced after
# library(frmtmb).
local({
  ns <- asNamespace("frmtmb")
  round0 <- function(frame, lp) {
    if (!is.null(lp[["constant"]]) || !is.null(lp[["equate"]]) ||
          !is.null(lp[["nl_body"]])) return(1)
    if (!identical(lp[["link"]]$name %||% "", "identity")) return(1)
    if (!identical(lp[["dpar"]], "mu")) return(1)
    fam <- frame[["spec"]]$responses[[lp[["resp"]]]]$family
    fname <- fam[["family"]]
    if (!is.character(fname) || length(fname) != 1L ||
          !fname %in% predictor_unit_families) return(1)
    init_fn <- fam[["init_dpars"]][["sigma"]]
    if (!is.function(init_fn)) return(1)
    sig <- Filter(function(l) {
      identical(l[["resp"]], lp[["resp"]]) && identical(l[["dpar"]], "sigma")
    }, frame[["linpreds"]])
    if (length(sig) != 1L || !is.null(sig[[1L]][["constant"]]) ||
          !is.null(sig[[1L]][["equate"]])) return(1)
    y <- frame[["y"]][[lp[["resp"]]]]
    if (!is.numeric(y) || !is.null(dim(y))) return(1)
    at <- frame[["aterm_values"]][[lp[["resp"]]]]
    w <- at[["weights"]]
    keep <- if (is.numeric(w) && length(w) == length(y) && any(w <= 0)) {
      which(w > 0)
    }
    raw <- tryCatch({
      res <- if (length(formals(init_fn)) >= 3L) {
        mu_residuals(frame, lp[["resp"]])
      }
      if (!is.null(keep)) {
        at <- subset_obs(at, keep, length(y))
        if (!is.null(res)) res <- res[keep]
        y <- y[keep]
      }
      if (is.null(res)) init_fn(y, at) else init_fn(y, at, res)
    }, error = function(e) NA_real_)
    raw <- suppressWarnings(as.numeric(raw)[1L])
    if (!is.finite(raw) || raw <= 0) return(1)
    2^round(log2(raw))
  }
  environment(round0) <- ns
  utils::assignInNamespace("predictor_unit_of", round0, "frmtmb")
  cat("predictor_unit_of: the round-0 rule (the start of sigma)\n")
})
