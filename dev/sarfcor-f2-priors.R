# Fix round 1, F2: brms's class priors on a registry term's bounded
# parameters (lagsar, errorsar), core side. Rscript <this>
p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/R/priors.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
# the class row in get_prior()/default_prior(), with its interval
rep1("  add <- function(class, coef = \"\", group = \"\", dpar = \"\", nlpar = \"\",
                  resp = \"\", lb = NA_real_, dflt = NULL) {",
"  add <- function(class, coef = \"\", group = \"\", dpar = \"\", nlpar = \"\",
                  resp = \"\", lb = NA_real_, dflt = NULL, ub = NA_real_) {")
rep1("      resp = resp, dpar = dpar, nlpar = nlpar, lb = lb,
      ub = NA_real_, source = \"default\"",
"      resp = resp, dpar = dpar, nlpar = nlpar, lb = lb,
      ub = ub, source = \"default\"")
rep1("    for (cl in names(autocor_prior_classes)) {
      if (length(autocor_class_idx(ac, cl))) add(cl, resp = resp_lab)
    }",
"    for (cl in names(autocor_prior_classes)) {
      if (length(autocor_class_idx(ac, cl))) add(cl, resp = resp_lab)
    }
    # a registry term's bounded parameters (sar()'s lagsar), under
    # their natural names and on their interval, as brms lists them
    for (pr in autocor_seam_bounded(ac)) {
      add(pr$name, resp = resp_lab, lb = pr$lower, ub = pr$upper)
    }")
# routing: a natural-scale spec that names a registry term's parameter
rep1("    bad_shape <- dpar_shape_refusal(fit, s)
    if (!is.null(bad_shape)) frm_stop(bad_shape, call. = FALSE)
    ord_th <- ",
"    if (isTRUE(s$natural) && s$class %in% autocor_seam_classes(frame)) {
      resolve_ac_class(s)
      next
    }
    bad_shape <- dpar_shape_refusal(fit, s)
    if (!is.null(bad_shape)) frm_stop(bad_shape, call. = FALSE)
    ord_th <- ")
rep1("autocor_class_idx <- function(ac, cls) {
  st <- ac[[\"struct\"]]",
"autocor_class_idx <- function(ac, cls) {
  if (autocor_is_seam(ac)) {
    hit <- Filter(function(pr) identical(pr$name, cls),
                  autocor_seam_bounded(ac))
    return(if (length(hit)) hit[[1L]]$idx)
  }
  st <- ac[[\"struct\"]]")
rep1("autocor_trans <- function(ac, cls) {",
"autocor_trans <- function(ac, cls) {
  if (autocor_is_seam(ac)) {
    hit <- Filter(function(pr) identical(pr$name, cls),
                  autocor_seam_bounded(ac))[[1L]]
    return(list(map = \"bounded\", lo = hit$lower, hi = hit$upper))
  }")
rep1("ac_trans_value <- function(th, tr) {
  if (identical(tr$map, \"identity\")) return(th)",
"ac_trans_value <- function(th, tr) {
  if (identical(tr$map, \"identity\")) return(th)
  if (identical(tr$map, \"bounded\")) {
    return(tr$lo + (tr$hi - tr$lo) / (1 + exp(-th[1])))
  }")
rep1("ac_trans_logjac <- function(th, tr) {
  if (identical(tr$map, \"identity\")) return(0)",
"ac_trans_logjac <- function(th, tr) {
  if (identical(tr$map, \"identity\")) return(0)
  if (identical(tr$map, \"bounded\")) {
    # d rho / d theta = (hi - lo) p (1 - p), the same form as cosy's
    s <- 1 / (1 + exp(-th[1]))
    return(log(tr$hi - tr$lo) + log(s) + log(1 - s))
  }")
rep1("ac_bound_theta <- function(v, tr, ac, cls, what) {
  if (identical(tr$map, \"identity\")) return(v)",
"ac_bound_theta <- function(v, tr, ac, cls, what) {
  if (identical(tr$map, \"identity\")) return(v)
  if (identical(tr$map, \"bounded\")) {
    if (v <= tr$lo || v >= tr$hi) {
      frm_stop(\"class = \\\"\", cls, \"\\\" \", what, \" = \", v, \" is outside \",
               \"the interval of \", ac[[\"label\"]], \", (\",
               format(tr$lo, digits = 6), \", \", format(tr$hi, digits = 6),
               \")\", call. = FALSE)
    }
    return(stats::qlogis((v - tr$lo) / (tr$hi - tr$lo)))
  }")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
cat("ok\n")
