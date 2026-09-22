# emmeans, insight and marginaleffects RUN on four fits, so that "no
# interop entry point breaks" is a measurement and not a grep.
#
# The same script runs against the base build and against the lane, and
# the two outputs are diffed:
#   Rscript dev/shapes-interop.R base > dev/shapes-log/interop-base.txt
#   Rscript dev/shapes-interop.R lane > dev/shapes-log/interop-lane.txt
#
# Punch round 1 rebuilt it. The first version printed
#   round(head(as.numeric(insight::get_residuals(f)), 3), 4)
# and as.numeric() flattens a matrix COLUMN FIRST, so when
# get_residuals() started returning the whole 150 x 4 summary matrix the
# first three numbers were still the first three residuals, bit for bit,
# and the diff was empty. A value probe cannot see a shape change. So
# every accessor now prints its CLASS, LENGTH, DIM and DIMNAMES before
# any value, and each one carries a contract that is CHECKED: the script
# prints FAIL lines and exits 1 when a contract breaks, so it is a
# harness and not a log.

arg <- commandArgs(trailingOnly = TRUE)
lib <- if (identical(arg[1], "base")) {
  "C:/Users/adf44/source/r/rellib-r3"
} else {
  "C:/Users/adf44/source/r/shapes-lib"
}
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("library:", lib, " frmtmb ",
    as.character(utils::packageVersion("frmtmb")), "\n")

nullish <- function(a, b) if (is.null(a)) b else a

FAILS <- new.env(parent = emptyenv())
FAILS$n <- 0L
FAILS$msg <- character(0)
check <- function(ok, what) {
  ok <- isTRUE(ok)
  cat("   [", if (ok) "ok  " else "FAIL", "] ", what, "\n", sep = "")
  if (!ok) {
    FAILS$n <- FAILS$n + 1L
    FAILS$msg <- c(FAILS$msg, what)
  }
  invisible(ok)
}

# class, length, dim and dimnames of whatever came back, printed before
# any value, because the shape is the thing a value probe cannot see
shape_line <- function(r) {
  cat("   class:  ", paste(class(r), collapse = ", "), "\n", sep = "")
  cat("   length: ", length(r), "\n", sep = "")
  cat("   dim:    ",
      if (is.null(dim(r))) "NULL" else paste(dim(r), collapse = " x "),
      "\n", sep = "")
  dn <- if (!is.null(dim(r))) dimnames(r) else list(names(r))
  txt <- vapply(nullish(dn, list(NULL)), function(v) {
    if (is.null(v)) "NULL" else {
      paste0(paste(utils::head(v, 6), collapse = ", "),
             if (length(v) > 6) " ...")
    }
  }, "")
  cat("   dimnames: ", paste(txt, collapse = " | "), "\n", sep = "")
  invisible(NULL)
}

show <- function(lab, expr, digits = 4, print_value = TRUE) {
  cat("\n-- ", lab, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(expr,
         warning = function(w) {
           cat("   [warning] ", conditionMessage(w), "\n", sep = "")
           invokeRestart("muffleWarning") },
         message = function(m) invokeRestart("muffleMessage")),
    error = function(e) {
      cat("   ERROR: ", conditionMessage(e), "\n", sep = "")
      structure(list(msg = conditionMessage(e)), class = "ioErr") })
  shape_line(r)
  if (inherits(r, "ioErr")) return(invisible(r))
  if (print_value) {
    v <- r
    if (is.numeric(v)) v <- round(v, digits)
    if (is.data.frame(v)) {
      for (j in seq_along(v)) {
        if (is.numeric(v[[j]])) v[[j]] <- round(v[[j]], digits)
      }
    }
    if (is.matrix(v) && nrow(v) > 6) v <- utils::head(v, 3)
    if (is.atomic(v) && is.null(dim(v)) && length(v) > 8) {
      v <- utils::head(v, 5)
    }
    print(v)
  }
  invisible(r)
}

set.seed(20260917)
n <- 150
dd <- data.frame(x = rnorm(n), z = rnorm(n),
                 f = factor(rep(c("a", "b", "c"), length.out = n)),
                 g = factor(rep(1:15, each = 10)))
dd$y <- rnorm(n, 1 + 0.5 * dd$x + c(a = 0, b = 0.4, c = -0.3)[dd$f], 1)
dd$bin <- rbinom(n, 1, plogis(0.2 + 0.7 * dd$x))
dd$ymix <- rnorm(n, 1 + 0.5 * dd$x + rnorm(15, 0, 0.6)[dd$g], 1)
dd$ord <- factor(cut(1 + 0.8 * dd$x + rnorm(n), c(-Inf, -0.3, 0.8, Inf),
                     labels = 1:3), ordered = TRUE)

fits <- list(
  gaussian = frm(bf(y ~ x + f) + gaussian(), data = dd),
  binomial = frm(bf(bin ~ x) + bernoulli(), data = dd),
  mixed = frm(bf(ymix ~ x + (1 | g)) + gaussian(), data = dd),
  ordinal = frm(bf(ord ~ x) + cumulative(), data = dd)
)

for (nm in names(fits)) {
  f <- fits[[nm]]
  cat("\n================ ", nm, " ================\n", sep = "")

  # emmeans
  show("emmeans(x)", {
    e <- emmeans::emmeans(f, "x", at = list(x = c(-1, 0, 1)))
    as.data.frame(summary(e))
  })
  if (nm == "gaussian") {
    show("emmeans pairwise f", {
      as.data.frame(summary(emmeans::emmeans(f, pairwise ~ f)$contrasts))
    })
  }

  # insight. Every accessor below is a VECTOR accessor in insight's own
  # contract: one number per observation. The length check is the thing
  # the value probe could not do.
  p <- show("insight::get_predicted", insight::get_predicted(f))
  check(!inherits(p, "ioErr") && length(as.numeric(p)) == n,
        "get_predicted() is one value per observation")
  pars <- show("insight::get_parameters",
               as.data.frame(insight::get_parameters(f)))
  show("insight::find_parameters", unlist(insight::find_parameters(f)))
  V <- show("insight::get_varcov diag", diag(insight::get_varcov(f)))
  check(!inherits(pars, "ioErr") && !inherits(V, "ioErr") &&
          length(V) == nrow(pars),
        "get_varcov() covers every get_parameters() row")
  show("insight::model_info family", insight::model_info(f)$family)
  rr <- show("insight::get_residuals", insight::get_residuals(f))
  check(!inherits(rr, "ioErr") && is.null(dim(rr)) && length(rr) == n,
        "get_residuals() is a length-n vector, not a summary matrix")
  check(inherits(rr, "insight_residuals"),
        "get_residuals() carries insight's class")
  rs <- show("insight::get_response", insight::get_response(f))
  check(!inherits(rs, "ioErr") && length(rs) == n,
        "get_response() is one value per observation")
  # fitted() and residuals() are brms's summary matrix on purpose, so
  # they carry no length contract; their shape is printed because the
  # insight accessors above are built on them and a change there is
  # what BLOCKER 1 was
  show("stats::fitted", stats::fitted(f))
  show("stats::residuals", stats::residuals(f))

  # the pairing lmtest::coeftest() makes: it keeps only the rows both
  # coef() and vcov() name, so a naming split drops rows in silence.
  # A fit WITH random effects has no such pairing in any of lme4,
  # glmmTMB or brms, since coef() is per-group there, so the contract
  # is asserted only where coef() is one vector.
  cf <- show("coef()", stats::coef(f))
  vn <- show("vcov() rownames", rownames(stats::vcov(f)))
  ct <- show("lmtest::coeftest", as.data.frame(unclass(lmtest::coeftest(f))))
  if (is.numeric(cf) && is.null(dim(cf))) {
    check(!inherits(ct, "ioErr") && !inherits(vn, "ioErr") &&
            nrow(ct) == length(vn),
          "coeftest() keeps every vcov() row")
    check(setequal(names(cf), vn),
          "coef() and vcov() name the same parameters")
  } else {
    cat("   [n/a ] coeftest(): coef() is per-group here\n")
  }

  # marginaleffects
  sl <- show("marginaleffects::avg_slopes", {
    s <- as.data.frame(marginaleffects::avg_slopes(f, variables = "x"))
    s[, intersect(c("term", "group", "estimate", "std.error"), names(s))]
  })
  check(!inherits(sl, "ioErr") && all(is.finite(sl$std.error)),
        "avg_slopes() standard errors are finite")
  show("marginaleffects::predictions head", {
    p <- as.data.frame(marginaleffects::predictions(
      f, newdata = marginaleffects::datagrid(x = c(-1, 0, 1))))
    p[, intersect(c("group", "x", "estimate", "std.error"), names(p))]
  })
}

# An OUTSIDE JUDGE for the one seam where the lane's own number was
# wrong and nothing in this script could tell: the ordinal average
# marginal effect. MASS::polr fits the same model with its own
# machinery, so its standard errors do not come through frmtmb.
cat("\n================  ordinal against MASS::polr  ================\n")
pol <- show("MASS::polr avg_slopes", {
  pm <- MASS::polr(ord ~ x, data = dd, Hess = TRUE, method = "logistic")
  s <- as.data.frame(marginaleffects::avg_slopes(pm, variables = "x"))
  s[, c("group", "estimate", "std.error")]
})
frm_sl <- show("frmtmb avg_slopes", {
  s <- as.data.frame(marginaleffects::avg_slopes(fits$ordinal,
                                                 variables = "x"))
  s[, c("group", "estimate", "std.error")]
})
if (!inherits(pol, "ioErr") && !inherits(frm_sl, "ioErr")) {
  m <- merge(pol, frm_sl, by = "group", suffixes = c(".polr", ".frm"))
  m$se_ratio <- round(m$std.error.frm / m$std.error.polr, 4)
  print(m)
  check(all(abs(m$se_ratio - 1) < 0.05),
        "ordinal avg_slopes() SE within 5% of MASS::polr on every category")
}

cat("\n---- contract summary ----\n")
cat("failures:", FAILS$n, "\n")
if (FAILS$n) {
  cat(paste0("  FAIL: ", FAILS$msg, collapse = "\n"), "\n")
  cat("\nHARNESS FAILED\n")
  quit(status = 1L)
}
cat("\nDONE\n")
