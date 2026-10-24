# Guardrail 1 probe. Usage: Rscript review-splitsurv-g1.R <arm>
# arm: base (rellib-r11 only), nosurv (lane lib, survival not attached),
# surv (lane lib, library(frmtmb.survival)), ns (namespace loaded only)
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm != "base") "C:/Users/adf44/source/r/wt-splitsurv-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
if (arm == "surv") suppressMessages(library(frmtmb.survival))
if (arm == "ns") loadNamespace("frmtmb.survival")
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("survival ns loaded:", isNamespaceLoaded("frmtmb.survival"), "\n")
set.seed(1)
dd <- data.frame(x = rnorm(150))
dd$time <- rexp(150, exp(-0.5 + 0.7 * dd$x))
dd$out <- rbinom(150, 1, 0.3)
try1 <- function(lab, expr) {
  r <- tryCatch({
    v <- expr
    if (inherits(v, "frmtmb_fit")) {
      paste0("FIT logLik=", format(as.numeric(logLik(v)), digits = 15))
    } else paste0("OK class=", paste(class(v), collapse = "/"))
  }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
  cat(sprintf("%-34s %s\n", lab, r))
}
f <- bf(time | cens(out) ~ x)
try1("frm(family = \"cox\")", frm(f, family = "cox", data = dd))
try1("frm(family = \"Cox\")", frm(f, family = "Cox", data = dd))
try1("brmsfamily(\"cox\")", brmsfamily("cox"))
try1("frm(brmsfamily(\"cox\"))", frm(f, family = brmsfamily("cox"), data = dd))
try1("frm(family = cox())", frm(f, family = cox(), data = dd))
try1("cox_baseline exists", get("cox_baseline"))
try1("frm(brms::cox())", frm(f, family = brms::cox(), data = dd))
try1("frm(brms::cox(bhaz=df6))",
     frm(f, family = brms::cox(bhaz = list(df = 6)), data = dd))
try1("frm(brms::brmsfamily(\"cox\"))",
     frm(f, family = brms::brmsfamily("cox"), data = dd))
try1("get_prior(family = \"cox\")",
     get_prior(f, data = dd, family = "cox"))
try1("frm(family = \"royston_parmar\")",
     frm(bf(time | cens(out) ~ x), family = "royston_parmar", data = dd))
try1("brmsfamily(\"royston_parmar\")", brmsfamily("royston_parmar"))
try1("frm(family = \"wiener\")", frm(f, family = "wiener", data = dd))
# frm_compat() rows for cox
cc <- tryCatch(frm_compat(), error = function(e) NULL)
if (!is.null(cc)) {
  cat("frm_compat rows:", nrow(cc), "\n")
}
