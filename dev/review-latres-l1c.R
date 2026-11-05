# Which step of the L1 check has no factor on the two cells, on the test's data?
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
Sys.setenv(FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/stan-cache", FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages(library(rstan))
hel <- new.env(parent = asNamespace("frmtmb"))
for (f in c("helper-brms.R", "helper-brms-methods.R", "helper-warnings.R")) sys.source(file.path("tests/testthat", f), envir = hel)
assignInNamespace("latres_refuse_margin", function(fit) invisible(NULL), ns = "frmtmb")
t0 <- parse(text = readLines("tests/testthat/test-latres-brms.R"))
for (e in t0) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% c("lrb_data")) eval(e)
d <- lrb_data()
box0 <- get("latres_disp_box", asNamespace("frmtmb"))
for (cc in list(list("Beta", brms::brmsfamily("beta"), Beta(), "unit"), list("Gamma()", brms::brmsfamily("Gamma"), Gamma(), "pos"))) {
  for (b in c(box0, Inf)) {
    assignInNamespace("latres_disp_box", b, ns = "frmtmb")
    f <- stats::as.formula(paste(cc[[4]], "~ x + unstr(t, g)"))
    fit <- suppressWarnings(frm(frmtmb::bf(f), family = cc[[3]], data = d))
    bk <- fi("frame_latres")(fit$frame)[[1]]
    th <- fit$estimates$theta[bk$theta_idx]
    C <- as.matrix(fi("us_chol_cor")(th[-1], bk$aux$d))
    r <- tryCatch(hel$brms_lp_check(brms::bf(f), cc[[2]], d, fit, joint = TRUE), error = function(e) conditionMessage(e))
    cat(cc[[1]], "box", b, " code", fit$opt$convergence, " ll", as.numeric(logLik(fit)), " log disp", fit$opt$par[names(fit$opt$par) == "betad"][1], " max theta", max(abs(th[-1])), " min eig", min(eigen(C, TRUE, TRUE)$values), "\n   check:", if (is.character(r)) substr(r, 1, 160) else paste("const", r$measured_const), "\n")
  }
}
