# Re-check 1: the mo() span argument behind the new branch of
# check_cs_identified(), and its false-alarm side.
#   Rscript dev/csfactor-rev2-mo.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== build", TAG, ":", as.character(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                           1L, 320L)
fitq <- function(...) {
  tryCatch(suppressWarnings(suppressMessages(frm(...))),
           error = function(e) structure(list(m = short(e)),
                                         class = "revfail"))
}
# one line per shape: refused with its message, or fitted with the df, the
# log likelihood, the count of non-finite standard errors and the smallest
# eigenvalue of vcov (negative or zero means a flat or worse direction)
rep1 <- function(tag, r) {
  if (inherits(r, "revfail")) {
    cat(sprintf("%-34s REFUSED  %s\n", tag, r$m)); return(invisible(NULL))
  }
  v <- tryCatch(suppressWarnings(vcov(r)), error = function(e) NULL)
  bad <- if (is.null(v)) NA_integer_ else sum(!is.finite(sqrt(diag(v))))
  ev <- if (is.null(v)) NA_real_ else
    min(eigen(v, symmetric = TRUE, only.values = TRUE)$values)
  cat(sprintf(
    "%-34s FITTED   df=%-3d logLik=%.9f  bad-se=%s  min-eig=%.4g\n",
    tag, attr(logLik(r), "df"), as.numeric(logLik(r)), bad, ev))
  invisible(r)
}

# --------------------------------------------------------------- data
# `m` gets a NON-UNIFORM monotone effect so the fitted simplex is far
# from uniform: the jump sits between the third and fourth level.
set.seed(3151)
n <- 600
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$zb <- rbinom(n, 1L, 0.5)
mcode <- sample(0:3, n, TRUE)
d$m <- factor(mcode, levels = 0:3, ordered = TRUE)
d$k <- mcode                                    # the integer spelling
d$kf <- factor(mcode)
# a DIFFERENT ordered factor, correlated with m but not a function of it
m2code <- pmin(3L, pmax(0L, mcode + sample(c(-1L, 0L, 1L), n, TRUE)))
d$m2 <- factor(m2code, levels = 0:3, ordered = TRUE)
# a COARSENING of m: two groups of two levels
d$mc <- factor(c("A", "A", "B", "B")[mcode + 1L])
# a REFINEMENT of m: each m level split in two
d$mr <- factor(paste0(mcode, sample(c("p", "q"), n, TRUE)))
# an unrelated factor
d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
d$g <- factor(sample(c("p", "q"), n, TRUE))
d$m1 <- d$m
eff <- c(0, 0.15, 0.30, 1.30)[mcode + 1L]
eta <- 0.5 * d$x + eff
p1 <- plogis(-0.7 - eta); p2 <- plogis(0.9 - eta)
P <- cbind(p1, p2 - p1, 1 - p2)
d$yo <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))
d$yo2 <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))
cat("table(m) =", paste(table(d$m), collapse = "/"),
    " table(mc) =", paste(table(d$mc), collapse = "/"),
    " nlevels(mr) =", nlevels(d$mr),
    " table(yo) =", paste(table(d$yo), collapse = "/"), "\n")
cat("cross-tab m against m2 (m2 is NOT a function of m):\n")
print(table(d$m, d$m2))

# ------------------------------------------- 1a the span, numerically
cat("\n### 1a  the fitted mo() column against the span of Dmo\n")
fmo <- fitq(bf(yo ~ mo(m)), family = sratio(), data = d)
if (inherits(fmo, "revfail")) {
  cat("mo(m) alone FAILED: ", fmo$m, "\n")
} else {
  lp <- fmo$frame$linpreds[["yo.mu"]]
  mi <- lp[["mo"]][[1L]]
  zr <- fmo$estimates[[mi[["zeta"]]]]
  zeta <- exp(c(0, zr)); zeta <- zeta / sum(zeta)
  cat("zeta = ", paste(sprintf("%.6f", zeta), collapse = " "),
      "   (uniform would be ", sprintf("%.6f", 1 / length(zeta)),
      " each)\n", sep = "")
  cat("max/min zeta ratio = ",
      sprintf("%.3f", max(zeta) / min(zeta)), "\n", sep = "")
  cz0 <- c(0, cumsum(zeta))
  bmo <- unname(fixef(fmo)["mom", "Estimate"])
  term <- bmo * mi[["D"]] * cz0[mi[["codes"]] + 1L]
  # the reference build has no such function; the basis is the same
  # definition either way, so fall back to it there
  Dm <- if (exists("mo_indicator_basis", asNamespace("frmtmb"))) {
    frmtmb:::mo_indicator_basis(mi)
  } else {
    B <- outer(as.integer(mi[["codes"]]), seq_len(mi[["D"]]), "==") * 1
    if (!is.null(mi[["mult"]])) B <- B * as.numeric(mi[["mult"]])
    colnames(B) <- paste0(mi[["label"]], ".", seq_len(mi[["D"]]))
    B
  }
  cat("Dmo columns: ", paste(colnames(Dm), collapse = ","),
      "   ncol = ", ncol(Dm), "\n", sep = "")
  res <- stats::lsfit(Dm, term, intercept = FALSE)$residuals
  cat("||term|| = ", sprintf("%.9g", sqrt(sum(term^2))),
      "   ||residual of term on Dmo|| = ", sprintf("%.6e", sqrt(sum(res^2))),
      "\n", sep = "")
  cat("term values at codes 0..3 = ",
      paste(sprintf("%.6f", tapply(term, mi[["codes"]], `[`, 1L)),
            collapse = " "), "\n", sep = "")
  # the span claim is about EVERY simplex, not this one: ten random
  # simplices, same projection
  set.seed(99)
  worst <- 0
  for (b in 1:10) {
    zz <- stats::rgamma(length(zeta), 0.4); zz <- zz / sum(zz)
    tt <- stats::rnorm(1) * mi[["D"]] * c(0, cumsum(zz))[mi[["codes"]] + 1L]
    rr <- stats::lsfit(Dm, tt, intercept = FALSE)$residuals
    worst <- max(worst, sqrt(sum(rr^2)))
  }
  cat("worst residual over 10 random simplices = ",
      sprintf("%.6e", worst), "\n", sep = "")
}

# -------------------------------------------- 1b the false-alarm side
cat("\n### 1b  refusals that should happen\n")
rep1("mo(m) + cs(m)", fitq(bf(yo ~ mo(m) + cs(m)), family = sratio(),
                           data = d))
rep1("mo(m) + cs(mr) REFINEMENT",
     fitq(bf(yo ~ mo(m) + cs(mr)), family = sratio(), data = d))
rep1("mo(k) integer + cs(kf)",
     fitq(bf(yo ~ mo(k) + cs(kf)), family = sratio(), data = d))
rep1("mo(m1) + mo(m2) + cs(m1)",
     fitq(bf(yo ~ mo(m1) + mo(m2) + cs(m1)), family = sratio(), data = d))
rep1("x + mo(m) + cs(m)",
     fitq(bf(yo ~ x + mo(m) + cs(m)), family = sratio(), data = d))

cat("\n### 1b  shapes that must NOT be refused\n")
rep1("mo(m) alone", fitq(bf(yo ~ mo(m)), family = sratio(), data = d))
rep1("cs(m) alone", fitq(bf(yo ~ cs(m)), family = sratio(), data = d))
rep1("mo(m) + cs(x)", fitq(bf(yo ~ mo(m) + cs(x)), family = sratio(),
                           data = d))
rep1("mo(m) + cs(f) unrelated",
     fitq(bf(yo ~ mo(m) + cs(f)), family = sratio(), data = d))
rep1("mo(m) + cs(mc) COARSENING",
     fitq(bf(yo ~ mo(m) + cs(mc)), family = sratio(), data = d))
rep1("mo(m) + cs(m2) correlated",
     fitq(bf(yo ~ mo(m) + cs(m2)), family = sratio(), data = d))
rep1("mo(m1) + mo(m2) + cs(x)",
     fitq(bf(yo ~ mo(m1) + mo(m2) + cs(x)), family = sratio(), data = d))
rep1("mo(m):z + cs(m)  z continuous",
     fitq(bf(yo ~ mo(m):z + cs(m)), family = sratio(), data = d))
rep1("mo(m):zb + cs(m)  z binary",
     fitq(bf(yo ~ mo(m):zb + cs(m)), family = sratio(), data = d))
rep1("mo(m)*z + cs(m)",
     fitq(bf(yo ~ mo(m) * z + cs(m)), family = sratio(), data = d))

cat("\n### 1b  is the accepted mo():z case really identified?\n")
for (tg in c("z", "zb")) {
  fo <- fitq(bf(as.formula(paste0("yo ~ mo(m):", tg, " + cs(m)"))),
             family = sratio(), data = d)
  fb <- fitq(bf(yo ~ cs(m)), family = sratio(), data = d)
  if (!inherits(fo, "revfail") && !inherits(fb, "revfail")) {
    cat("  mult = ", tg, ": logLik ", sprintf("%.7f", as.numeric(logLik(fo))),
        " at df ", attr(logLik(fo), "df"), " against cs(m) alone ",
        sprintf("%.7f", as.numeric(logLik(fb))), " at df ",
        attr(logLik(fb), "df"), "\n", sep = "")
    cat("  gain = ", sprintf("%.6f", as.numeric(logLik(fo)) -
                               as.numeric(logLik(fb))),
        " for ", attr(logLik(fo), "df") - attr(logLik(fb), "df"),
        " df\n", sep = "")
    v <- suppressWarnings(vcov(fo))
    ev <- eigen(v, symmetric = TRUE, only.values = TRUE)$values
    cat("  vcov eigenvalues: ", paste(sprintf("%.4g", ev), collapse = " "),
        "\n", sep = "")
    cat("  all standard errors finite: ", all(is.finite(sqrt(diag(v)))),
        "   vcov positive definite: ", all(ev > 0), "\n", sep = "")
    print(signif(fixef(fo), 4))
  }
}

cat("\n### 1b  thres(gr = ) and multivariate: which refusal wins\n")
r <- fitq(bf(yo | thres(gr = g) ~ mo(m) + cs(m)), family = sratio(),
          data = d)
cat("thres(gr) + mo(m) + cs(m): ",
    if (inherits(r, "revfail")) r$m else "FITTED", "\n", sep = "")
r <- fitq(bf(yo ~ mo(m) + cs(m)) + bf(yo2 ~ x + cs(f)),
          family = sratio(), data = d)
cat("mv, bad lp first: ",
    if (inherits(r, "revfail")) r$m else "FITTED", "\n", sep = "")
r <- fitq(bf(yo ~ x + cs(f)) + bf(yo2 ~ mo(m) + cs(m)),
          family = sratio(), data = d)
cat("mv, bad lp second: ",
    if (inherits(r, "revfail")) r$m else "FITTED", "\n", sep = "")

cat("\n### 1b  the hole the SECOND condition leaves open\n")
# mo(fo) beside fo itself: X already spans the codes, so the mo branch
# declines to fire by design. Is that model identified?
rep1("mo(m) + m + cs(x)",
     fitq(bf(yo ~ mo(m) + m + cs(x)), family = sratio(), data = d))
rep1("mo(m) + m (no cs at all)",
     fitq(bf(yo ~ mo(m) + m), family = sratio(), data = d))

cat("\n### 1b  the refusal with a proper prior, for consistency\n")
r <- fitq(bf(yo ~ mo(m) + cs(m)), family = sratio(), data = d,
          prior = set_prior("normal(0, 1)", class = "b"))
cat("mo(m) + cs(m) with a normal(0,1) prior on b: ",
    if (inherits(r, "revfail")) "REFUSED" else "FITTED", "\n", sep = "")
cat("\nDONE ", TAG, "\n")
