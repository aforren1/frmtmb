## Reviewer, claim 1: the SAME model in real brms (dev/brmsnames-rev-brms.R)
## and in frmtmb.sample (dev/brmsnames-rev-frm.R), accessor by accessor:
## class, names, dims and dimnames, never values. No shim.
##   Rscript dev/brmsnames-rev-shapes.R [C1 ...]
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample)); q(library(brms))
av <- commandArgs(trailingOnly = TRUE)
cases <- if (length(av)) av else paste0("C", 1:8)

shape <- function(x, depth = 0) {
  if (inherits(x, "error")) return(paste("ERROR"))
  d <- list(class = setdiff(class(x), c("frmtmb_hypothesis")))
  if (!is.null(dim(x))) d$dim_len <- length(dim(x))
  if (!is.null(dimnames(x))) {
    dn <- dimnames(x)
    d$dimnames <- lapply(dn, function(v) if (is.null(v)) NULL else v)
    d$dimnames_names <- names(dn)
    if (length(dim(x)) >= 2) d$dim_tail <- dim(x)[-1]
  }
  if (is.list(x) && !is.data.frame(x) && depth < 3) {
    d$names <- names(x)
    d$el <- lapply(x, shape, depth = depth + 1)
  }
  if (is.data.frame(x)) { d$names <- names(x) }
  d
}
diffshape <- function(a, b, path = "") {
  out <- character(0)
  if (is.character(a) && is.character(b) && length(a) == 1 && length(b) == 1 &&
      (a == "ERROR" || b == "ERROR")) {
    if (!identical(a, b)) out <- c(out, paste0(path, ": brms ", if (is.list(a)) "ok" else a,
                                              ", frmtmb ", if (is.list(b)) "ok" else b))
    return(out)
  }
  if (!is.list(a) || !is.list(b)) {
    if (!identical(a, b)) {
      fa <- paste(utils::head(setdiff(a, b), 6), collapse = " ")
      fb <- paste(utils::head(setdiff(b, a), 6), collapse = " ")
      if (!nzchar(fa) && !nzchar(fb)) {
        out <- c(out, paste0(path, ": same set, different ORDER or dup (brms ",
                             length(a), ", frmtmb ", length(b), ")"))
      } else {
        out <- c(out, paste0(path, ": brms-only [", fa, "] frmtmb-only [", fb, "]"))
      }
    }
    return(out)
  }
  for (k in union(names(a), names(b))) {
    if (is.null(k) || !nzchar(k)) next
    out <- c(out, diffshape(a[[k]], b[[k]], paste0(path, "/", k)))
  }
  if (is.null(names(a)) && length(a) && length(a) == length(b)) {
    for (i in seq_along(a)) out <- c(out, diffshape(a[[i]], b[[i]], paste0(path, "[", i, "]")))
  }
  out
}

t1 <- function(e) tryCatch(q(e), error = function(err) err)
calls <- function(x, is_brms, nm) {
  grp <- "g"
  eff <- if (nm == "C5") "az" else if (nm == "C6") "sigma_z" else "x"
  list(
    variables = t1(sort(grep("^(b_|sd_|cor_|r_|sigma$)", variables(x), value = TRUE))),
    variables_dups = t1({v <- variables(x); unique(v[duplicated(v)])}),
    fixef = t1(fixef(x)),
    fixef_nosum = t1(fixef(x, summary = FALSE)),
    ranef = t1(ranef(x)),
    ranef_nosum = t1(ranef(x, summary = FALSE)),
    coef = t1(coef(x)),
    coef_nosum = t1(coef(x, summary = FALSE)),
    VarCorr = t1(VarCorr(x)),
    VarCorr_nosum = t1(VarCorr(x, summary = FALSE)),
    posterior_summary_cols = t1(colnames(posterior_summary(x))),
    posterior_summary_b = t1(posterior_summary(x, variable = "^b_", regex = TRUE)),
    summary_cols = t1({s <- summary(x); if (is_brms) colnames(s$fixed) else colnames(s)}),
    as_draws_df_cols = t1(sort(intersect(grep("^(b_|r_)", names(as_draws_df(x)), value = TRUE),
                                         grep("^(b_|r_)", names(as_draws_df(x)), value = TRUE)))),
    as_array = t1(as.array(x, variable = "b_Intercept")),
    as_matrix = t1(as.matrix(x, variable = c("b_Intercept"))),
    hyp = t1(hypothesis(x, "Intercept > 0")),
    hyp_ranef = t1(hypothesis(x, "Intercept > 0", scope = "ranef", group = grp)),
    hyp_coef = t1(hypothesis(x, "Intercept > 0", scope = "coef", group = grp)),
    bayes_R2 = t1(bayes_R2(x)),
    log_lik = t1(dim(log_lik(x, ndraws = 5))[2]),
    ce = t1(names(conditional_effects(x, effects = eff)[[1]]))
  )
}

for (nm in cases) {
  bp <- sprintf("dev/stan-cache/brmsnames-rev-brms-%s.rds", nm)
  fp <- sprintf("dev/stan-cache/brmsnames-rev-frm-%s.rds", nm)
  if (!file.exists(bp) || !file.exists(fp)) { cat("==", nm, "missing input\n"); next }
  bf_ <- readRDS(bp); fo <- readRDS(fp)
  cat("\n==", nm, "==\n")
  if (inherits(bf_, "error")) { cat("brms refused the model:", conditionMessage(bf_), "\n")
    cat("frmtmb variables(fit):", variables(fo$fit), "\n"); next }
  if (inherits(fo$ds, "error")) { cat("frmtmb draws missing\n"); next }
  a <- calls(bf_, TRUE, nm); b <- calls(fo$ds, FALSE, nm)
  for (k in names(a)) {
    sa <- shape(a[[k]]); sb <- shape(b[[k]])
    if (inherits(a[[k]], "error") || inherits(b[[k]], "error")) {
      ea <- if (inherits(a[[k]], "error")) substr(conditionMessage(a[[k]]), 1, 90) else "ok"
      eb <- if (inherits(b[[k]], "error")) substr(conditionMessage(b[[k]]), 1, 90) else "ok"
      if (!identical(ea, eb)) cat(sprintf("  %-22s brms: %s | frmtmb: %s\n", k, ea, eb))
      next
    }
    if (k %in% c("variables", "posterior_summary_cols", "summary_cols", "as_draws_df_cols",
                 "variables_dups", "log_lik", "ce")) {
      if (!identical(a[[k]], b[[k]])) {
        cat(sprintf("  %-22s brms-only: %s\n  %-22s frmtmb-only: %s\n", k,
                    paste(utils::head(setdiff(a[[k]], b[[k]]), 12), collapse = " "), "",
                    paste(utils::head(setdiff(b[[k]], a[[k]]), 12), collapse = " ")))
        if (!length(setdiff(a[[k]], b[[k]])) && !length(setdiff(b[[k]], a[[k]])))
          cat("     (same set; order or duplicates differ)\n")
      }
      next
    }
    ds <- diffshape(sa, sb, k)
    if (length(ds)) cat(paste0("  ", utils::head(ds, 8), collapse = "\n"), "\n")
  }
}
cat("DONE\n")
