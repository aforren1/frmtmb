## REVIEW claim 2: the pin is applied ONLY when the subset's remaining
## category labels are an INITIAL SEGMENT of the fitted ones. A subset
## that loses an INTERIOR category is left unpinned. What does
## influence() then return for that row: a misaligned coefficient vector,
## an NA, or a refusal?
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, " ver =",
    as.character(utils::packageVersion("frmtmb")), "\n\n")
options(width = 140)

## n = 60, four categories, EXACTLY ONE row in category 2 (interior).
mk <- function(seed = 901, n = 60) {
  set.seed(seed)
  x <- stats::rnorm(n)
  cp <- cbind(stats::plogis(-0.7 - 0.5 * x), stats::plogis(0.6 - 0.5 * x),
              stats::plogis(2.0 - 0.5 * x))
  y <- 1L + rowSums(stats::runif(n) > cp)
  mid <- which(y == 2L)
  if (length(mid) < 2L) stop("seed gives too few category-2 rows")
  ## move every category-2 row but one down to category 1
  y[mid[-1L]] <- 1L
  if (!any(y == 4L)) stop("seed reached no top category")
  data.frame(x = x, y = y)
}
dd <- mk()
imid <- which(dd$y == 2L)
cat("table(y) =", paste(table(dd$y), collapse = "/"),
    " single interior row =", imid, "\n\n")

show <- function(tag, inf, imid) {
  cat("--", tag, "\n")
  cat("   cols:", paste(colnames(inf$fixed), collapse = ","), "\n")
  cat("   NA cells in whole table:", sum(is.na(inf$fixed)), "\n")
  cat("   NA cells in the interior row:", sum(is.na(inf$fixed[imid, ])), "\n")
  cat("   interior row:",
      paste(signif(inf$fixed[imid, ], 9), collapse = "  "), "\n")
  cd <- suppressWarnings(cooks.distance(inf))
  cat("   cooks.distance at that row:", cd[imid],
      "  NA in cooks:", sum(is.na(cd)), "\n")
}

for (code in c("integer", "ordered", "unordered", "character")) {
  d <- dd
  d$y <- switch(code,
                integer = dd$y,
                ordered = factor(dd$y, levels = 1:4, ordered = TRUE),
                unordered = factor(dd$y, levels = 1:4),
                character = as.character(dd$y))
  fit <- try(suppressWarnings(frm(bf(y ~ x), family = cumulative(),
                                  data = d)), silent = TRUE)
  if (inherits(fit, "try-error")) {
    cat("==", code, ": frm() REFUSED:",
        sub("\n.*", "", conditionMessage(attr(fit, "condition"))), "\n\n")
    next
  }
  cat("==", code, ": n_tau =", length(fit$estimates[["tau_raw"]]), "\n")
  inf <- try(suppressWarnings(influence(fit, force = TRUE)), silent = TRUE)
  if (inherits(inf, "try-error")) {
    cat("   influence() REFUSED:",
        conditionMessage(attr(inf, "condition")), "\n\n")
    next
  }
  show(code, inf, imid)
  ## what the honest answer looks like: the same model pinned by hand on
  ## the subset with the interior category recoded back to 1..4
  sub <- dd[-imid, , drop = FALSE]
  ref <- try(suppressWarnings(frm(bf(y | thres(3) ~ x),
                                  family = cumulative(), data = sub)),
             silent = TRUE)
  if (!inherits(ref, "try-error")) {
    w <- frmtmb:::get_coef.frmtmb_fit(ref)
    cat("   thres(3) hand reference (integer codes 1,3,4 kept as is):\n      ",
        paste(names(w), collapse = "  "), "\n      ",
        paste(signif(w, 9), collapse = "  "), "\n")
    cat("   logLik pinned reference =", sprintf("%.10f", logLik(ref)), "\n")
  } else {
    cat("   thres(3) hand reference REFUSED:",
        sub("\n.*", "", conditionMessage(attr(ref, "condition"))), "\n")
  }
  ## and the recount a bare frm() on the subset gives
  bare <- try(suppressWarnings(frm(bf(y ~ x), family = cumulative(),
                                   data = sub)), silent = TRUE)
  if (!inherits(bare, "try-error")) {
    w <- frmtmb:::get_coef.frmtmb_fit(bare)
    cat("   bare frm() on the subset: n_tau =",
        length(bare$estimates[["tau_raw"]]), " coefs:",
        paste(paste0(names(w), "=", signif(w, 9)), collapse = "  "), "\n")
    cat("   logLik bare =", sprintf("%.10f", logLik(bare)), "\n")
  }
  cat("\n")
}

## ---- grouped, interior loss: is the row SHIFTED? -------------------
cat("#### grouped thres(gr = g) with a single interior-category row\n")
mkg <- function(seed = 902, n = 96) {
  set.seed(seed)
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- stats::rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.0), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
  y <- vapply(seq_len(n), function(i) {
    1L + sum(stats::runif(1) >
               stats::plogis(tau[[as.character(g[i])]] - 0.5 * x[i]))
  }, 1L)
  mid <- which(y == 2L)
  y[mid[-1L]] <- 1L
  data.frame(x = x, g = g, y = y)
}
dg <- mkg()
img <- which(dg$y == 2L)
cat("table(y) =", paste(table(dg$y), collapse = "/"),
    " per-level max =", paste(tapply(dg$y, dg$g, max), collapse = "/"),
    " interior row =", img, " in level", as.character(dg$g[img]), "\n")
for (code in c("integer", "ordered")) {
  d <- dg
  if (identical(code, "ordered")) {
    d$y <- factor(dg$y, levels = 1:4, ordered = TRUE)
  }
  fit <- try(suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                                  family = cumulative(), data = d)),
             silent = TRUE)
  if (inherits(fit, "try-error")) {
    cat("==", code, ": frm() REFUSED:",
        sub("\n.*", "", conditionMessage(attr(fit, "condition"))), "\n")
    next
  }
  th <- fit$spec$responses$y$family[["thres"]]
  cat("==", code, ": nthres =", paste(th[["nthres"]], collapse = "/"),
      " nraw =", length(fit$estimates[["tau_raw"]]), "\n")
  inf <- suppressWarnings(influence(fit, force = TRUE))
  show(code, inf, img)
  ## the honest pinned reference
  sub <- dg[-img, , drop = FALSE]
  sub$k <- as.integer(th[["nthres"]][match(as.character(sub$g),
                                           th[["levels"]])])
  ref <- try(suppressWarnings(frm(bf(y | thres(k, gr = g) ~ x),
                                  family = cumulative(), data = sub)),
             silent = TRUE)
  if (!inherits(ref, "try-error")) {
    w <- frmtmb:::get_coef.frmtmb_fit(ref)
    cat("   thres(k, gr) reference:\n      ",
        paste(signif(w, 9), collapse = "  "), "\n")
  } else {
    cat("   reference REFUSED:",
        sub("\n.*", "", conditionMessage(attr(ref, "condition"))), "\n")
  }
  cat("\n")
}
cat("DONE rev-02\n")
