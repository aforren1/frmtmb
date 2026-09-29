# Every consumer of a cs() factor term, on the fixed build: fitting,
# fitted(), frm_linpred(), predict(), simulate(), conditional_effects(),
# emmeans(), summary()/fixef()/variables(), default_prior(),
# set_prior(class = "b"), frm_sample() names. Seed 405.
#   Rscript dev/csfactor-consumers.R <lib> > dev/csfactor-log/consumers.txt
args <- commandArgs(TRUE)
LIB <- if (length(args)) args[[1L]] else
  "C:/Users/adf44/source/r/wt-csfactor-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "\n")

set.seed(405)
n <- 500
x <- rnorm(n)
fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(fc)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n)
d <- data.frame(x, fc, yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
d$fb <- as.numeric(d$fc == "b")
d$fcc <- as.numeric(d$fc == "c")

ff <- frm(bf(yo ~ cs(fc)), family = sratio(), data = d)
fd <- frm(bf(yo ~ cs(fb) + cs(fcc)), family = sratio(), data = d)
nd1 <- data.frame(fc = factor("c"))
nd1d <- data.frame(fb = 0, fcc = 1)
nd3 <- data.frame(fc = factor(c("a", "b", "c")))

go <- function(lab, expr) {
  cat("\n--", lab, "--\n")
  print(tryCatch(expr, error = function(e) paste("ERROR:",
                                                 conditionMessage(e))))
}

go("frm_linpred response, 3 rows",
   round(frm_linpred(ff, newdata = nd3, type = "response"), 6))
go("frm_linpred response, 1 row factor('c')",
   round(frm_linpred(ff, newdata = nd1, type = "response"), 6))
go("hand dummies, same row",
   round(frm_linpred(fd, newdata = nd1d, type = "response"), 6))
go("fitted() se at 1 row",
   round(fitted(ff, newdata = nd1), 6))
go("fitted() in sample, max |diff| against hand dummies",
   max(abs(fitted(ff)[, "Estimate", ] - fitted(fd)[, "Estimate", ])))
go("residuals() in sample, max |diff|",
   max(abs(residuals(ff)[, "Estimate"] - residuals(fd)[, "Estimate"])))
go("predict() at factor('c'), 20000 draws", {
  set.seed(3)
  round(predict(ff, newdata = nd1, ndraws = 20000L,
                propagate_error = FALSE), 6)
})
go("simulate() at factor('c'), 5000 draws", {
  set.seed(4)
  prop.table(table(simulate(ff, newdata = nd1, nsim = 5000L)))
})
go("conditional_effects(effects = 'fc')", {
  ce <- conditional_effects(ff, effects = "fc")
  head(as.data.frame(ce[[1L]])[, c("fc", "cats__", "estimate__")], 9)
})
go("emmeans(~ fc)", {
  if (!requireNamespace("emmeans", quietly = TRUE)) "emmeans absent" else {
    summary(emmeans::emmeans(ff, ~ fc))
  }
})
go("summary() fixed block", summary(ff)$fixed)
go("fixef()", round(fixef(ff), 6))
go("variables()", variables(ff))
go("default_prior()", default_prior(ff))
go("hypothesis() on a bcs row",
   hypothesis(ff, "bcs_fcc[1] = 0")$hypothesis)
go("set_prior class b coef fcc shrinks bcs_fcc", {
  fp <- frm(bf(yo ~ cs(fc)), family = sratio(), data = d,
            prior = set_prior("normal(0, 0.05)", class = "b",
                              coef = "fcc"))
  rbind(flat = c(ff$estimates$bcs3, ff$estimates$bcs4),
        shrunk = c(fp$estimates$bcs3, fp$estimates$bcs4))
})
go("set_prior class b (no coef) covers all four", {
  fp <- frm(bf(yo ~ cs(fc)), family = sratio(), data = d,
            prior = set_prior("normal(0, 0.05)", class = "b"))
  unlist(fp$estimates[c("bcs1", "bcs2", "bcs3", "bcs4")])
})
go("frm_sample() names", {
  if (!requireNamespace("frmtmb.sample", quietly = TRUE)) {
    "frmtmb.sample absent"
  } else {
    dr <- frmtmb.sample::frm_sample(bf(yo ~ cs(fc)), family = sratio(),
                                    data = d, chains = 1L, iter = 400L,
                                    refresh = 0L, seed = 7L)
    grep("^bcs", frmtmb::variables(dr), value = TRUE)
  }
})
go("new level at newdata is refused",
   frm_linpred(ff, newdata = data.frame(fc = factor("zz")),
               type = "response"))
go("control: no false alarm on z + cs(x)", {
  d$z <- rnorm(n)
  as.numeric(logLik(frm(bf(yo ~ z + cs(x)), family = sratio(), data = d)))
})
go("refusal: fc + cs(fc)",
   frm(bf(yo ~ fc + cs(fc)), family = sratio(), data = d))
cat("\ndone\n")
