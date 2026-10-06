# Defect 8 on a given build: per seed of brms_monotonic's data code,
# `ls ~ mo(income) * age` and `ls ~ mo(income)`: the fit-time warnings,
# how many standard errors are finite, which are lost and why, whether
# conditional_effects(fit, "income:age") prints, and the coefficient
# SEs against a reference that holds every saturated simplex coordinate
# (weight below 1e-8) fixed: the inverse of the exact Hessian over the
# other parameters.
#   Rscript dev/nanse-mo-lane.R [lib] [seeds] [out]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/wt-nanse-lib"
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:200
out <- if (length(args) > 2) args[3] else "dev/nanse-log/mo-lane.tsv"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
grDevices::pdf(NULL)
mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
catch <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = v, warnings = w)
}
fmt <- function(x) paste(signif(x, 4), collapse = ",")
rows <- list()
for (s in seeds) {
  d <- mk(s)
  for (form in c("int", "main")) {
    fml <- if (form == "int") ls ~ mo(income) * age else ls ~ mo(income)
    r <- catch(frm(fml, data = d))
    f <- r$value
    s2 <- catch(sdr <- frmtmb:::sdr_of(f))
    nm <- frmtmb:::outer_par_names(f)
    se <- suppressWarnings(sqrt(diag(sdr$cov.fixed)))
    b <- nm %in% c("(Intercept)", "age", "moincome", "moincome:age",
                   "sigma_(Intercept)")
    # the reference: hold the saturated simplex coordinates
    p <- f$opt$par
    H <- f$obj$he(p)
    zt <- f$estimates[grepl("^zeta", names(f$estimates))]
    sat <- rep(FALSE, length(p))
    zpos <- which(grepl("^zeta", nm))
    w_all <- unlist(lapply(zt, function(z) {
      x <- exp(c(0, z))
      x / sum(x)
    }))
    for (zn in names(zt)) {
      x <- exp(c(0, zt[[zn]]))
      x <- x / sum(x)
      pos <- which(startsWith(nm, paste0(zn, "_")))
      # a weight k > 1 at 0 is its own coordinate; the reference weight
      # at 0 moves every coordinate of the block together
      sat[pos[x[-1] < 1e-8]] <- TRUE
      if (x[1] < 1e-8) sat[pos] <- TRUE
    }
    keep <- !sat
    Vr <- tryCatch(solve(H[keep, keep]), error = function(e) NULL)
    se_ref <- rep(NA_real_, length(p))
    if (!is.null(Vr)) se_ref[keep] <- suppressWarnings(sqrt(diag(Vr)))
    ce <- if (form == "int") {
      tryCatch({
        x <- suppressWarnings(conditional_effects(f, "income:age"))
        print(x)
        "OK"
      }, error = function(e) conditionMessage(e))
    } else NA_character_
    lost <- sdr$se_lost
    rows[[length(rows) + 1L]] <- data.frame(
      seed = s, form = form, code = f$opt$convergence,
      n_warn_fit = length(r$warnings),
      se_warn = any(grepl("Standard errors are not available", r$warnings)),
      warn_fit = substr(paste(r$warnings, collapse = " | "), 1, 300),
      n_warn_sdr = length(s2$warnings),
      se_finite = sum(is.finite(se)), n_par = length(se),
      n_lost = length(lost),
      lost = paste(paste0(names(lost), "=", lost), collapse = ","),
      n_sat = sum(w_all < 1e-8),
      ce_ok = identical(ce, "OK") || is.na(ce),
      ce_msg = if (is.na(ce)) "" else substr(ce, 1, 80),
      se_b = fmt(se[b]), se_b_ref = fmt(se_ref[b]),
      max_rel_b = max(abs(se[b] / se_ref[b] - 1)),
      stringsAsFactors = FALSE)
  }
  if (s %% 20 == 0) cat("seed", s, "\n")
}
X <- do.call(rbind, rows)
utils::write.table(X, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", nrow(X), "rows to", out, "\n")
