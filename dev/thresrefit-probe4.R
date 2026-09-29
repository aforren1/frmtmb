## Probe 4 on the REFERENCE build: what breaks downstream of an
## imputation set whose threshold counts differ, and whether B is
## reachable through frmtmb.sample's prior-predictive path.
lib <- Sys.getenv("FRMTMB_PROBE_LIB",
                  "C:/Users/adf44/source/r/rellib-r3")
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("frmtmb", as.character(packageVersion("frmtmb")), "from", lib, "\n")
say <- function(...) cat(..., "\n", sep = "")
try_msg <- function(expr) {
  tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
}

mk <- function(seed, n = 60, tau = c(-0.6, 0.5, 2.6), slope = 0.5) {
  set.seed(seed)
  x <- rnorm(n)
  cp <- cbind(plogis(tau[1] - slope * x), plogis(tau[2] - slope * x),
              plogis(tau[3] - slope * x))
  u <- runif(n)
  data.frame(x = x, z = rnorm(n), y = 1L + rowSums(u > cp))
}
dd <- mk(404)
d2 <- dd
d2$y[d2$y == 4L] <- 3L
say("A5: imputation tables: ",
    paste(table(factor(dd$y, 1:4)), collapse = "/"), " and ",
    paste(table(factor(d2$y, 1:4)), collapse = "/"))
imps <- list(dd, d2, dd)
m <- suppressWarnings(frm_multiple(bf(y ~ x + z), data = imps,
                                   family = cumulative()))
say("A5: per-fit tau_raw lengths = ",
    paste(vapply(m$fits, function(f) length(f$estimates$tau_raw), 1L),
          collapse = ","))
say("A5: pooled rows = ", paste(rownames(m$pooled), collapse = ","))
s <- try_msg(summary(m))
say("A5: summary() class = ", paste(class(s), collapse = "/"))
if (is.character(s) && length(s) == 1L && grepl("^ERROR", s)) say(s)

m0 <- suppressWarnings(frm_multiple(bf(y ~ x), data = imps,
                                    family = cumulative()))
for (meth in c("D1", "D2", "D3")) {
  r <- try_msg(anova(m, m0, method = meth))
  if (is.character(r)) say("A5: anova(", meth, ") -> ", r) else {
    say("A5: anova(", meth, ") ok:")
    print(r)
  }
}
h <- try_msg(hypothesis(m, "x = 0"))
say("A5: hypothesis -> ",
    if (is.character(h)) h else "ok")

## same imputations where every one keeps all four categories: the
## control that says the breakage is about the count, not the pooling
imps_ok <- list(dd, dd, dd)
mo <- suppressWarnings(frm_multiple(bf(y ~ x + z), data = imps_ok,
                                    family = cumulative()))
mo0 <- suppressWarnings(frm_multiple(bf(y ~ x), data = imps_ok,
                                     family = cumulative()))
for (meth in c("D1", "D2", "D3")) {
  r <- try_msg(anova(mo, mo0, method = meth))
  say("A5 control: anova(", meth, ") -> ",
      if (is.character(r)) r else "ok")
}

## ------------------------------------------------------------------
## B through frmtmb.sample: the prior-predictive path
## ------------------------------------------------------------------
ok <- requireNamespace("frmtmb.sample", quietly = TRUE)
say("B: frmtmb.sample available = ", ok)
if (ok) {
  fns <- ls(asNamespace("frmtmb.sample"))
  say("B: frmtmb.sample names mentioning prior: ",
      paste(grep("prior", fns, value = TRUE), collapse = ", "))
}
say("B: core names calling draw_prior_pars:")
for (nm in ls(asNamespace("frmtmb"), all.names = TRUE)) {
  f <- get(nm, envir = asNamespace("frmtmb"))
  if (!is.function(f)) next
  txt <- paste(deparse(body(f)), collapse = " ")
  if (grepl("draw_prior_pars", txt, fixed = TRUE)) say("  ", nm)
}
say("B: names calling prior_entry_label:")
for (nm in ls(asNamespace("frmtmb"), all.names = TRUE)) {
  f <- get(nm, envir = asNamespace("frmtmb"))
  if (!is.function(f)) next
  txt <- paste(deparse(body(f)), collapse = " ")
  if (grepl("prior_entry_label", txt, fixed = TRUE)) say("  ", nm)
}
