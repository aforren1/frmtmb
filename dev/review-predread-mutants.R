# Item 1 and 9: mutant builds of core, one per defect a reader could carry.
#   Rscript dev/review-predread-mutants.R <name> : source copy under
#   rev-predread-src/<name>, patched once (a pattern must match exactly once),
#   installed into rev-predread-lib/<name>.
src <- "C:/Users/adf44/source/r/frmtmb-wt-predread"
SRC <- "C:/Users/adf44/source/r/rev-predread-src"
LIB <- "C:/Users/adf44/source/r/rev-predread-lib"
muts <- list(
  cs_se = list("R/predict.R",
    "for (ct in ord_cs_values(object, lp, newdata, n)) {\n      ep <- which(jc$names == ct$par)",
    "for (ct in list()) {\n      ep <- which(jc$names == ct$par)"),
  cs_se_rev = list("R/predict.R",
    "B[(k - 1L) * n + seq_len(n), k] <- ct$vals",
    "B[(k - 1L) * n + seq_len(n), K1 - k + 1L] <- ct$vals"),
  cs_factor_nd = list("R/predict.R",
    "v <- if (is.null(cols)) ct[[\"vals\"]] else cols[[i]]",
    "v <- if (is.null(cols)) ct[[\"vals\"]] else if (all(cols[[i]] %in% c(0, 1))) cols[[i]] * 0 else cols[[i]]"),
  me_var_double = list("R/predict.R",
    "for (mp in ed[[\"me_parts\"]]) out <- out + (mp$g * mp$sd)^2",
    "for (mp in ed[[\"me_parts\"]]) out <- out + 2 * (mp$g * mp$sd)^2"),
  me_var_drop = list("R/predict.R",
    "for (mp in ed[[\"me_parts\"]]) out <- out + (mp$g * mp$sd)^2",
    "for (mp in ed[[\"me_parts\"]]) out <- out + 0"),
  me_shared_z = list("R/me.R",
    "drawn <- stats::rnorm(length(first), o$x[first], o$sd[first])",
    "drawn <- o$x[first] + o$sd[first] * stats::rnorm(1)"),
  me_half_sd = list("R/me.R",
    "drawn <- stats::rnorm(length(first), o$x[first], o$sd[first])",
    "drawn <- stats::rnorm(length(first), o$x[first], o$sd[first] / 2)"),
  me_gr_perrow = list("R/me.R",
    "J <- match(gc, lev)\n    me_check_level_values(xv, sv, J, lev, t)",
    "J <- match(gc, lev)\n    me_check_level_values(xv, sv, J, lev, t)\n    J <- seq_along(xv)"),
  link_identity = list("R/predict.R",
    "if (!identical(fam[[\"type\"]], \"ordinal\")) return(lp[[\"link\"]])",
    "return(lp[[\"link\"]])"),
  emm_me_design = list("R/interop.R",
    "!length(lp[[\"gps\"]]) && !length(lp[[\"mo\"]]) && !length(lp[[\"mi\"]]) &&\n    !length(lp[[\"me\"]])",
    "!length(lp[[\"gps\"]]) && !length(lp[[\"mo\"]]) && !length(lp[[\"mi\"]])"),
  me_store_off = list("R/me.R",
    "if (!is.null(old) && identical(old$x, xv)) return(old$value)",
    "if (FALSE) return(old$value)"),
  curve_check_ratio = list("R/curve-cov.R",
    "rel <- max(abs(se[!miss] - se_ref[!miss]) /\n               pmax(se_ref[!miss], .Machine$double.eps))",
    "rel <- max(abs(se[!miss] / pmax(se_ref[!miss], .Machine$double.eps) - 1))"),
  resp_read_off = list("R/predict.R",
    "resp_read <- function(object, resp) {\n  if (length(object$spec$responses) == 1L) NULL else resp\n}",
    "resp_read <- function(object, resp) {\n  resp\n}"),
  me_share = list("R/me.R",
    "key <- paste(lp[[\"resp\"]], lp[[\"dpar\"]], k, sep = \"\\r\")",
    "key <- as.character(k)"),
  band_var = list("R/me.R",
    "extra <- ifelse(gr > 0 & mm$s2 > 0, (as.numeric(se_coef) / gr)^2 / mm$s2, 0)",
    "extra <- ifelse(mm$var > 0, as.numeric(se_coef)^2 / mm$var, 0)"),
  band_cat = list("R/predict.R",
    "  at <- if (!is.null(mm)) matrix(mm$plug, dim(out)[1L], dim(out)[3L])\n",
    paste0(
      "  if (!is.null(mm)) {\n",
      "    zs <- stats::qnorm(stats::ppoints(99L))\n",
      "    vals <- lapply(zs, mm$val)\n",
      "    vr <- matrix(mm$var, dim(out)[1L], dim(out)[3L])\n",
      "    for (k in seq_len(dim(out)[3L])) {\n",
      "      V <- matrix(vapply(vals, function(v) as.matrix(v)[, k], numeric(dim(out)[1L])), dim(out)[1L], length(zs))\n",
      "      m <- as.numeric(out[, 1L, k]) / up\n",
      "      se <- as.numeric(out[, 2L, k]) / up\n",
      "      infl <- sqrt(1 + pmax(se^2 - vr[, k] / up^2, 0) / pmax(vr[, k] / up^2, .Machine$double.eps))\n",
      "      for (i in seq_along(q)) {\n",
      "        v <- apply(V, 1L, stats::quantile, probs = probs[i], na.rm = TRUE, names = FALSE) / up\n",
      "        out[, q[i], k] <- stats::plogis(stats::qlogis(m) + (stats::qlogis(v) - stats::qlogis(m)) * infl) * up\n",
      "      }\n",
      "    }\n",
      "    return(out)\n",
      "  }\n",
      "  at <- NULL\n")),
  emm_cs_design = list("R/interop.R",
    "design <- !epred && !use_re && !any(cs_K1 > 0L) &&",
    "design <- !epred && !use_re &&")
)
nm <- commandArgs(TRUE)[1]
if (nm == "list") {
  cat(names(muts), sep = "\n")
  quit("no")
}
m <- muts[[nm]]
stopifnot(!is.null(m))
d <- file.path(SRC, nm)
unlink(d, recursive = TRUE)
dir.create(d, recursive = TRUE)
for (p in c("DESCRIPTION", "NAMESPACE", ".Rbuildignore")) {
  file.copy(file.path(src, p), file.path(d, p))
}
for (p in c("R", "inst", "data", "man")) {
  file.copy(file.path(src, p), d, recursive = TRUE)
}
f <- file.path(d, m[[1L]])
x <- paste(readLines(f, warn = FALSE), collapse = "\n")
parts <- strsplit(x, m[[2L]], fixed = TRUE)[[1L]]
if (length(parts) != 2L) stop("pattern matched ", length(parts) - 1L, " times: ", nm)
writeLines(paste0(parts[1L], m[[3L]], parts[2L]), f, useBytes = TRUE)
lib <- file.path(LIB, nm)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_LIBS = paste("C:/Users/adf44/source/r/rellib-r16",
                          "C:/Users/adf44/AppData/Local/R/win-library/4.6", sep = ";"))
r <- system2(file.path(R.home("bin"), "R.exe"),
             c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "--no-test-load",
               paste0("--library=", lib), shQuote(d)), stdout = TRUE, stderr = TRUE)
cat(tail(r, 2), sep = "\n")
unlink(d, recursive = TRUE)
