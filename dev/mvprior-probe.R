# Lane wt-mvprior: what frmtmb does with each specification in
# dev/mvprior-cases.R, on the arm MVPRIOR_ARM names (base = rellib-r3,
# lane = the lane library). For each one: the internal parameters the
# resolver puts a density or a bound on, or the error and its class.
#
# Usage: MVPRIOR_ARM=base Rscript dev/mvprior-probe.R \
#          > dev/mvprior-log/probe-base.txt
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
suppressMessages(library(frmtmb))
source(file.path(MVPRIOR_ROOT, "dev/mvprior-cases.R"))
cat("arm", mvprior_arm, " frmtmb", format(packageVersion("frmtmb")),
    "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")

frm_spec <- function(s) {
  set_prior(s$prior, class = s$class, coef = s$coef, group = s$group,
            resp = s$resp, dpar = s$dpar, nlpar = s$nlpar, lb = s$lb)
}
lab <- function(s) {
  f <- c(class = s$class, coef = s$coef, group = s$group, resp = s$resp,
         dpar = s$dpar, nlpar = s$nlpar)
  f <- f[nzchar(f)]
  paste0(s$prior, if (!is.na(s$lb)) paste0(" lb=", s$lb), " [",
         paste(names(f), f, sep = "=", collapse = ", "), "]")
}
flat <- function(x) gsub("[[:space:]]+", " ", x)

reached <- function(design, pl) {
  r <- ns$resolve_priorlist(design, pl)
  pt <- design$frame[["par_template"]]
  nm <- vapply(r$entries, function(e) {
    nms <- ns$par_template_names(pt[[e$comp]], e$comp)
    paste(nms[e$idx], collapse = "+")
  }, "")
  c(if (length(nm)) paste(nm, collapse = " ") else "(no density)",
    if (length(r$lower)) paste0("lb:", paste(names(r$lower), collapse = ",")))
}

for (case in names(mvprior_cases)) {
  cc <- mvprior_cases[[case]]
  m <- mvprior_model(case, "frmtmb")
  cat("\n== ", case, ": ", cc$what, " ==\n", sep = "")
  design <- tryCatch(ns$prior_design(m[[1]], m[[2]], m[[3]], list()),
                     error = function(e) {
                       cat("MODEL REFUSED:", flat(conditionMessage(e)), "\n")
                       NULL
                     })
  if (is.null(design)) next
  for (s in cc$specs) {
    r <- tryCatch({
      pl <- frm_spec(s)
      paste("ACCEPT:", paste(reached(design, pl), collapse = " || "))
    }, error = function(e) {
      paste0("REFUSE [", paste(class(e)[1:2], collapse = ","), "]: ",
             flat(conditionMessage(e)))
    })
    cat(sprintf("%-58s %s\n", lab(s), r))
  }
  cat("-- default_prior() --\n")
  dp <- as.data.frame(default_prior(m[[1]], data = m[[2]], family = m[[3]]))
  dp <- dp[dp$class != "theta", ]
  print(dp[, c("prior", "class", "coef", "group", "resp", "dpar",
               "nlpar")], row.names = FALSE)
}
