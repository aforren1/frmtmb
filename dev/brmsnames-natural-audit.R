## Punch round 2, item 1: which distributional parameters can take the
## natural-scale flag of brms_coef_table(), and what brms does with each.
## The flag needs a non-primary dpar (it lands in `betad`) that is
## intercept-only with no formula written, so every non-primary dpar of
## every family can reach it. Listed per family: the dpar, its link, how
## the family reads it (elementwise through the link, or jointly as a
## multinomial logit), and brms's class for the same family and name.
##   Rscript dev/brmsnames-natural-audit.R > dev/brmsnames-log/natural-audit.txt
source("dev/brmsnames-libs.R")
brmsnames_libs("lane")
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
q(library(brms))
ns <- asNamespace("frmtmb")
link_name <- function(l) {
  if (is.character(l)) return(l)
  if (is.list(l) && !is.null(l$name)) return(l$name)
  "closure"
}
brms_dpars <- function(fname) {
  f <- tryCatch(q(brms::brmsfamily(fname)), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  tryCatch(brms:::valid_dpars(f), error = function(e) NULL)
}
rows <- list()
add <- function(pkg, label, fam, bname = NULL) {
  prim <- fam$primary_dpars %||% "mu"
  resp <- fam$post$dpar_response$dpars %||% character(0)
  bd <- if (!is.null(bname)) brms_dpars(bname)
  for (dp in setdiff(fam$dpars, prim)) {
    mix <- dp %in% resp && grepl("^theta[0-9]+$", dp)
    ls <- dp %in% fam$link_scale_dpars
    rows[[length(rows) + 1L]] <<- data.frame(
      pkg = pkg, family = label, dpar = dp,
      link = link_name(fam$links[[dp]]),
      read = if (mix) "softmax, all K" else if (ls) "softmax, one row" else
        "elementwise",
      unmodeled_is = if (mix) "simplex theta1..K" else if (ls)
        "coefficient b_" else "natural, linkinv",
      brms = if (identical(bname, "mixture")) {
        if (mix) "same class (simplex)" else "same class"
      } else
        if (is.null(bd)) "no brms family" else
        if (dp %in% bd) "same class" else "not a brms dpar",
      stringsAsFactors = FALSE)
  }
}
reg <- get("family_registry", ns)
for (nm in unique(names(reg))) {
  fam <- tryCatch(q(reg[[nm]]()), error = function(e) NULL)
  if (is.null(fam)) {
    cat("core family", nm, "needs arguments; built below if listed\n")
    next
  }
  add("frmtmb", nm, fam, nm)
}
add("frmtmb", "mixture(gaussian, gaussian)",
    frmtmb::mixture(gaussian(), gaussian()), "mixture")
add("frmtmb", "mixture(gaussian, gaussian, gaussian)",
    frmtmb::mixture(gaussian(), gaussian(), gaussian()), "mixture")
mv <- tryCatch(q(frmtmb::mixture_mvn(2, 2)), error = function(e) e)
if (inherits(mv, "frmtmb_family")) add("frmtmb", "mixture_mvn(2, 2)", mv)
if (inherits(mv, "error")) cat("mixture_mvn:", conditionMessage(mv), "\n")
ext <- list(
  frmtmb.eam = list(wiener = quote(wiener()), lba = quote(lba(2)),
                    rdm = quote(rdm(2)), gddm = quote(gddm()),
                    wiener_gng = quote(wiener_gng())),
  frmtmb.latent = list(`hmm(2)` = quote(hmm(2)), `lca(2)` = quote(lca(2))),
  frmtmb.learn = list(bandit2arm_delta = quote(bandit2arm_delta()),
                      bandit2arm_dual = quote(bandit2arm_dual()),
                      prl_fictitious = quote(prl_fictitious()),
                      bandit4arm2_kalman_filter =
                        quote(bandit4arm2_kalman_filter()),
                      ts_par7 = quote(ts_par7()),
                      igt_pvl_delta = quote(igt_pvl_delta()),
                      igt_orl = quote(igt_orl()), rlddm = quote(rlddm())))
for (pkg in names(ext)) {
  q(library(pkg, character.only = TRUE))
  for (lab in names(ext[[pkg]])) {
    fam <- tryCatch(q(eval(ext[[pkg]][[lab]], asNamespace(pkg))),
                    error = function(e) e)
    if (inherits(fam, "error")) {
      cat(pkg, lab, "not built:", conditionMessage(fam), "\n")
      next
    }
    add(pkg, lab, fam, if (lab == "wiener") "wiener")
  }
}
tab <- do.call(rbind, rows)
options(width = 200)
print(tab, row.names = FALSE)
cat("\njoint rows:\n")
print(tab[tab$read != "elementwise", ], row.names = FALSE)
cat("DONE\n")
