.libPaths(c(Sys.getenv("R2_LIB", "C:/Users/adf44/source/r/predfix-lib"),
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

# capture the pre-fit's mapped-back template, to measure whether the map
# back is exact: an exact map lands the template at the unscaled optimum
.r2 <- new.env()
suppressMessages(trace("autoscale_prefit", where = asNamespace("frmtmb"),
      exit = quote(assign("tpl", returnValue(), envir = .r2)),
      print = FALSE))

# fit and collect warnings/errors
fitw <- function(...) {
  .r2$tpl <- NULL
  w <- character(0)
  f <- tryCatch(withCallingHandlers(frm(...), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(fit = f, warn = w, tpl = .r2$tpl)
}

rel <- function(a, b) {
  a <- unlist(a); b <- unlist(b)
  if (length(a) != length(b)) return(NA_real_)
  max(abs(a - b) / pmax(abs(b), 1e-8))
}
absd <- function(a, b) {
  a <- unlist(a); b <- unlist(b)
  if (length(a) != length(b)) return(NA_real_)
  max(abs(a - b))
}

vc_num <- function(f) {
  v <- tryCatch(VarCorr(f), error = function(e) NULL)
  if (is.null(v)) return(NA)
  unlist(lapply(v, function(x) {
    if (is.list(x)) unlist(lapply(x, function(y) if (is.numeric(y)) y))
    else if (is.numeric(x)) x
  }))
}

summ <- function(label, a, b) {
  fa <- a$fit; fb <- b$fit
  if (inherits(fa, "error") || inherits(fb, "error")) {
    cat(sprintf("%-40s ERROR a=%s | b=%s\n", label,
                if (inherits(fa, "error")) conditionMessage(fa) else "ok",
                if (inherits(fb, "error")) conditionMessage(fb) else "ok"))
    return(invisible())
  }
  ll <- c(as.numeric(logLik(fa)), as.numeric(logLik(fb)))
  re_a <- tryCatch(unlist(ranef(fa)), error = function(e) NA)
  re_b <- tryCatch(unlist(ranef(fb)), error = function(e) NA)
  va <- tryCatch(vcov(fa), error = function(e) NA)
  vb <- tryCatch(vcov(fb), error = function(e) NA)
  pa <- tryCatch(predict(fa), error = function(e) NA)
  pb <- tryCatch(predict(fb), error = function(e) NA)
  cat(sprintf(paste0("%-40s engaged=%s/%s dLL=%.3e ll=%.9f fixef=%.2e ",
                     "ranef=%.2e VarCorr=%.2e vcov=%.2e se=%.2e pred=%.2e ",
                     "code=%s/%s nwarn=%d/%d\n"),
              label, !is.null(a$tpl), !is.null(b$tpl), ll[1] - ll[2], ll[2],
              rel(fixef(fa), fixef(fb)), rel(re_a, re_b),
              rel(vc_num(fa), vc_num(fb)), rel(va, vb),
              rel(sqrt(diag(va)), sqrt(diag(vb))),
              rel(pa, pb), fa$opt$convergence, fb$opt$convergence,
              length(a$warn), length(b$warn)))
}

# how far the final fit moved from the mapped-back template, per component
tpl_move <- function(r) {
  if (is.null(r$tpl) || inherits(r$fit, "error")) return("no template")
  est <- r$fit$estimates
  out <- vapply(intersect(names(r$tpl), c("beta", "betad", "theta", "b")),
                function(nm) absd(r$tpl[[nm]], est[[nm]]), 0)
  paste(names(out), sprintf("%.2e", out), collapse = " ")
}
