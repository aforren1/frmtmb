# Reviewer probe: family = "wiener", wiener() and brms::wiener() under
# different attach orders. Arg 1: "eam_last", "brms_last", "ns_only".
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mode <- commandArgs(TRUE)[1]
suppressMessages({
  library(frmtmb)
  if (mode == "eam_last") { library(brms); library(frmtmb.eam) }
  if (mode == "brms_last") { library(frmtmb.eam); library(brms) }
  if (mode == "ns_only") loadNamespace("frmtmb.eam")
})
cat("mode", mode, "; search():", paste(search()[2:5], collapse = " "), "\n")
set.seed(20261013)
d <- frmtmb.eam:::ddm_simulate(200, mu = 0.8, bs = 1.3, ndt = 0.3)
show <- function(lab, fam) {
  r <- tryCatch({
    fit <- suppressMessages(frm(frmtmb::bf(rt | dec(upper) ~ 1), family = fam, data = d))
    f <- family(fit)
    sprintf("logLik %.8f; links: %s", as.numeric(logLik(fit)),
            paste(names(f$links), vapply(f$links, function(l) l$name, ""),
                  sep = "=", collapse = " "))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-34s %s\n", lab, r))
}
if (mode != "ns_only") {
  cat("environment of wiener():", environmentName(environment(wiener)), "\n")
  show("family = wiener()", wiener())
}
show("family = \"wiener\"", "wiener")
show("family = brms::wiener()", brms::wiener())
show("brms::wiener(link_bs = 'identity')", brms::wiener(link_bs = "identity"))
show("frmtmb.eam::wiener()", frmtmb.eam::wiener())
show("frmtmb.eam::wiener(link_ndt='log')", frmtmb.eam::wiener(link_ndt = "log"))
show("brmsfamily(\"wiener\")", brmsfamily("wiener"))
r <- tryCatch(frm(frmtmb::bf(rt | dec(upper) ~ 1, ndt ~ 1), family = brms::wiener(),
                  data = transform(d, s = gl(4, 50))),
              error = function(e) conditionMessage(e))
cat("brms::wiener() with ndt ~ 1:", if (is.character(r)) r else
      sprintf("logLik %.8f", as.numeric(logLik(r))), "\n")
