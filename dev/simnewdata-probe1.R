# Probe: is brms fixture 1 (arma residual correlation) a structured
# draw, and what do the pieces do at newdata on it?
#   Rscript dev/simnewdata-probe1.R
source("dev/simnewdata-prelude.R")
suppressMessages(library(frmtmb))
d <- as.data.frame(brms:::brmsfit_example1$data)
attr(d, "terms") <- NULL
fit <- suppressWarnings(suppressMessages(frm(
  bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) +
       (1 + Trt | visit) + arma(visit, patient, cov = TRUE), sigma ~ Trt),
  family = student(), data = d)))
rspec <- frmtmb:::single_response(fit, "x")
cat("structured:", sim_is_structured(sim_context(fit, rspec, list(),
                                                 aterms = list())), "\n")
cat("names(fit):", names(fit), "\n")
cat("has $data:", !is.null(fit$data), "\n")
nd <- fit$data[1:10, ]
print(nd[, c("visit", "patient", "Trt", "Age")])
for (dn in names(rspec$dpars)) {
  r <- tryCatch(frm_linpred(fit, newdata = nd, dpar = dn,
                            type = "response"),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(dn, ":", format(r, digits = 4), "\n")
}
cat("blocks:", vapply(fit$frame$re_blocks, `[[`, "", "covstruct"), "\n")
