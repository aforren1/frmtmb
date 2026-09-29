# Reviewer, claim 2: which blocks re_formula keeps, on constructions the
# worker's probe does not build.
#
# The check is INDEPENDENT of lp_eta_design(): the expected population
# eta is rebuilt here as
#     eta(NULL) - sum over NON-smooth blocks of Z[, c_idx] %*% cvec[c_idx]
# straight off the frame, and compared with frm_linpred(re_formula = NA).
# The kept block ids are reported as well, for the record.
#   REVLIB=base Rscript dev/resmooth-rev-probe2.R
ARM <- if (identical(Sys.getenv("REVLIB"), "base")) "base" else "lane"
LIB <- if (ARM == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM:", ARM, "| frmtmb from:", find.package("frmtmb"), "\n\n")
ns <- asNamespace("frmtmb")
gt <- function(nm) get(nm, envir = ns)
msg <- function(e) substr(conditionMessage(e), 1, 150)
tryv <- function(x) tryCatch(x, error = function(e) e)

set.seed(5)
n <- 300
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                g = factor(rep(1:10, each = 30)),
                h = factor(rep(1:6, length.out = n)))
d$y <- sin(2 * pi * d$x) + d$z^2 + c(0, 1, -1)[d$f] * d$x +
  rnorm(10, 0, 0.5)[d$g] + rnorm(n, 0, 0.3)
d$y2 <- 0.7 * d$x + rnorm(n, 0, 0.4)
d$ys <- sin(2 * pi * d$x) + rnorm(n, 0, exp(-1 + 1.5 * d$x))
d$cnt <- rpois(n, exp(0.3 + sin(2 * pi * d$x)))

# kept block ids, as lp_eta_design reports them
kept_ids <- function(fit, use_re, nd = NULL) {
  bl <- fit$frame[["re_blocks"]]
  b1 <- vapply(bl, function(b) b[["b_idx"]][1], 0)
  k <- integer(0)
  for (lp in fit$frame[["linpreds"]]) {
    ed <- tryv(gt("lp_eta_design")(fit, lp, nd, use_re, FALSE))
    if (inherits(ed, "error")) return(paste0("ERR:", msg(ed)))
    k <- c(k, match(vapply(ed[["sm_blocks"]] %||% list(),
                           function(b) b[["b_idx"]][1], 0), b1))
    if (!is.null(nd) && length(ed[["sm_parts"]])) {
      k <- c(k, match(vapply(ed[["sm_parts"]],
                             function(p) p$bk[["b_idx"]][1], 0), b1))
    }
    if (use_re) {
      for (rp in ed[["re_parts"]] %||% list()) {
        k <- c(k, match(rp$bk[["b_idx"]][1], b1))
      }
    }
  }
  paste(sort(unique(k)), collapse = ",")
}

# the independent reconstruction of the population eta
pop_eta_ref <- function(fit, resp = NULL, dpar = "mu") {
  fr <- fit$frame
  cvec <- gt("expand_b")(fr, fit$estimates[["b"]], fit$estimates[["theta"]])
  lp <- NULL
  for (l in fr[["linpreds"]]) {
    if (!identical(l[["dpar"]], dpar)) next
    if (!is.null(resp) && !identical(l[["resp"]], resp)) next
    lp <- l; break
  }
  if (is.null(lp) || is.null(lp[["Z"]])) return(NULL)
  full <- as.vector(frm_linpred(fit, re_formula = NULL, resp = resp,
                                dpar = dpar, type = "link"))
  drop <- numeric(length(full))
  for (bk in fr[["re_blocks"]]) {
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp")) next
    ci <- bk[["c_idx"]]
    contrib <- as.numeric(lp[["Z"]][, ci, drop = FALSE] %*% cvec[ci])
    if (length(contrib) == length(drop)) drop <- drop + contrib
  }
  full - drop
}

probe <- function(label, form, family = gaussian(), resp = NULL,
                  dpar = "mu", gvar = "g") {
  fit <- tryv(suppressWarnings(frm(form, family = family, data = d)))
  if (inherits(fit, "error")) {
    cat(sprintf("%-26s FIT ERROR: %s\n", label, msg(fit))); return(NULL)
  }
  bl <- fit$frame[["re_blocks"]]
  cs <- vapply(bl, function(b) b[["covstruct"]], "")
  cat(sprintf("%-26s blocks %s\n", label,
              paste(seq_along(cs), cs, sep = ":", collapse = ",")))
  cat(sprintf("   kept: NA in-sample %-9s NA newdata %-9s NULL %-9s\n",
              kept_ids(fit, FALSE), kept_ids(fit, FALSE, d),
              kept_ids(fit, TRUE)))
  ena <- tryv(as.vector(frm_linpred(fit, re_formula = NA, resp = resp,
                                    dpar = dpar, type = "link")))
  enl <- tryv(as.vector(frm_linpred(fit, re_formula = NULL, resp = resp,
                                    dpar = dpar, type = "link")))
  ref <- tryv(pop_eta_ref(fit, resp, dpar))
  bar <- stats::as.formula(paste0("~ (1 | ", gvar, ")"))
  epf <- tryv(as.vector(frm_linpred(fit, re_formula = bar, resp = resp,
                                    dpar = dpar, type = "link")))
  cat(sprintf("   max|NA-NULL| %9.3g | max|NA - independent ref| %9.3g\n",
              if (inherits(ena, "error") || inherits(enl, "error")) NA_real_
              else max(abs(ena - enl)),
              if (inherits(ena, "error") || inherits(ref, "error") ||
                    is.null(ref)) NA_real_ else max(abs(ena - ref))))
  cat(sprintf("   re_formula = ~(1|%s): %s\n", gvar,
              if (inherits(epf, "error")) paste("ERROR:", msg(epf)) else
                sprintf("OK max|-NULL| %.3g max|-NA| %.3g",
                        max(abs(epf - enl)), max(abs(epf - ena)))))
  invisible(fit)
}

cat("---- 1. s(g, bs = 're') beside (1 | g) on the same factor ----\n")
probe("s(x)+s(g,re)+(1|g)", bf(y ~ s(x) + s(g, bs = "re") + (1 | g)))

cat("\n---- 2. an fs smooth with a by = factor ----\n")
probe("s(x,g,fs,by=f)", bf(y ~ s(x, g, bs = "fs", k = 5, by = f)))

cat("\n---- 3. t2(x, z, g, c(cr,cr,re)) ----\n")
probe("t2(x,z,g,cr,cr,re)", bf(y ~ t2(x, z, g, bs = c("cr", "cr", "re"))))

cat("\n---- 4. s(x, by = g): a by-FACTOR smooth, not group-indexed ----\n")
probe("s(x, by = g)", bf(y ~ s(x, by = g)))

cat("\n---- 5. gp(x, by = g) ----\n")
probe("gp(x, by = g)", bf(y ~ gp(x, by = g)))

cat("\n---- 6. te(x, z) ----\n")
probe("te(x, z)", bf(y ~ te(x, z)))

cat("\n---- 7. a smooth in a dpar formula: sigma ~ s(g, bs = 're') ----\n")
probe("sigma ~ s(g, re)", bf(ys ~ s(x), sigma ~ s(g, bs = "re")),
      dpar = "sigma")

cat("\n---- 8. multivariate, a smooth on ONE response ----\n")
fmv <- tryv(suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5)) +
                                   bf(y2 ~ x + (1 | g)),
                                 family = gaussian(), data = d)))
if (inherits(fmv, "error")) {
  cat("mv FIT ERROR:", msg(fmv), "\n")
} else {
  cs <- vapply(fmv$frame[["re_blocks"]], function(b) b[["covstruct"]], "")
  cat("mv blocks:", paste(seq_along(cs), cs, sep = ":", collapse = ","), "\n")
  for (rp in c("y", "y2")) {
    a <- as.vector(frm_linpred(fmv, re_formula = NA, resp = rp,
                               type = "link"))
    b <- as.vector(frm_linpred(fmv, re_formula = NULL, resp = rp,
                               type = "link"))
    r <- pop_eta_ref(fmv, rp, "mu")
    cat(sprintf("  resp %-3s max|NA-NULL| %9.3g | max|NA-ref| %9.3g\n", rp,
                max(abs(a - b)), if (is.null(r)) NA_real_ else
                  max(abs(a - r))))
  }
}

cat("\n---- 9. a mixture with a smooth ----\n")
fmix <- tryv(suppressWarnings(
  frm(bf(y ~ s(x, g, bs = "fs", k = 5)),
      family = mixture(gaussian(), gaussian()), data = d)))
if (inherits(fmix, "error")) {
  cat("mixture FIT ERROR:", msg(fmix), "\n")
} else {
  cs <- vapply(fmix$frame[["re_blocks"]], function(b) b[["covstruct"]], "")
  cat("mixture blocks:", paste(seq_along(cs), cs, sep = ":", collapse = ","),
      "\n")
  for (dp in c("mu1", "mu2")) {
    a <- tryv(as.vector(frm_linpred(fmix, re_formula = NA, dpar = dp,
                                    type = "link")))
    b <- tryv(as.vector(frm_linpred(fmix, re_formula = NULL, dpar = dp,
                                    type = "link")))
    r <- tryv(pop_eta_ref(fmix, NULL, dp))
    cat(sprintf("  dpar %-4s max|NA-NULL| %9.3g | max|NA-ref| %9.3g\n", dp,
                if (inherits(a, "error")) NA_real_ else max(abs(a - b)),
                if (inherits(a, "error") || inherits(r, "error") ||
                      is.null(r)) NA_real_ else max(abs(a - r))))
  }
}

cat("\n---- 10. the partial-formula refusal, three shapes ----\n")
fp <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | f) +
                                (1 | h)), family = gaussian(), data = d))
for (rf in list(~ (1 | f), ~ (1 | f) + (1 | h), ~ 0, ~ 1)) {
  r <- tryv(as.vector(frm_linpred(fp, re_formula = rf, type = "link")))
  cat(sprintf("  re_formula = %-22s %s\n", deparse1(rf),
              if (inherits(r, "error")) paste("ERROR:", msg(r)) else
                sprintf("OK max|-NULL| %.4g",
                        max(abs(r - as.vector(frm_linpred(fp,
                          re_formula = NULL, type = "link")))))))
}
fp1 <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | f)),
                            family = gaussian(), data = d))
r <- tryv(as.vector(frm_linpred(fp1, re_formula = ~ (1 | f), type = "link")))
cat("  ONE bar term, re_formula = ~(1|f):",
    if (inherits(r, "error")) paste("ERROR:", msg(r)) else
      sprintf("OK max|-NULL| %.4g",
              max(abs(r - as.vector(frm_linpred(fp1, re_formula = NULL,
                                                type = "link"))))), "\n")
cat("\nDONE\n")
