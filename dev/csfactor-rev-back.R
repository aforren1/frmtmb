# Two loose ends.
# 1. The backward-compatibility fallback: a frame with no cs_mm, which is
#    what a 0.64.0 object holds. Constructed by FITTING on the reference
#    build, saving, and predicting on the lane build.
# 2. Is there ANY multivariate cs() path that fails, which would give the
#    ?frm sentence a kernel of truth?
#   Rscript dev/csfactor-rev-back.R <lib> <tag> <rdsfile> <mode>
#     mode "save": fit and saveRDS; mode "load": readRDS and predict
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]; RDS <- args[[3L]]; MODE <- args[[4L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== build", TAG, ":", as.character(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), " mode", MODE, "\n")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                           1L, 250L)
go <- function(lab, expr) {
  cat("\n-- ", lab, " --\n", sep = "")
  print(tryCatch(expr, error = function(e) paste("ERROR:", short(e))))
}
set.seed(88)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(d$fc)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n); u2 <- runif(n)
d$yo <- ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L))
d$yo2 <- ifelse(u2 < p1, 1L, ifelse(u2 < p1 + p2, 2L, 3L))
nd <- data.frame(x = c(-1, 0, 1), z = 0,
                 fc = factor(c("a", "b", "c"), levels = c("a","b","c")))

if (identical(MODE, "save")) {
  ff <- frm(bf(yo ~ cs(x)), family = sratio(), data = d)
  cat("cs entries: ", length(ff$frame$linpreds[["yo.mu"]][["cs"]]),
      "  cs_mm: ", length(ff$frame$linpreds[["yo.mu"]][["cs_mm"]] %||%
                            list()), "\n", sep = "")
  go("frm_linpred on the SAVING build",
     round(frm_linpred(ff, newdata = nd, type = "response"), 6))
  saveRDS(ff, RDS)
  cat("saved\n")
} else {
  ff <- readRDS(RDS)
  lp <- ff$frame$linpreds[["yo.mu"]]
  cat("loaded. cs entries: ", length(lp[["cs"]]),
      "  cs_mm entries: ", length(lp[["cs_mm"]] %||% list()),
      "  first cs names: ", paste(names(lp[["cs"]][[1L]]), collapse = ","),
      "\n", sep = "")
  go("frm_linpred on the LANE build from a cs_mm-less frame",
     round(frm_linpred(ff, newdata = nd, type = "response"), 6))
  go("fitted() in sample from the same object",
     round(utils::head(fitted(ff)[, "Estimate", ], 3L), 6))
}

cat("\n### multivariate cs(): every prediction path\n")
mv <- suppressWarnings(suppressMessages(
  frm(bf(yo ~ x + cs(fc)) + bf(yo2 ~ z + cs(fc)), family = sratio(),
      data = d)))
mvn <- suppressWarnings(suppressMessages(
  frm(bf(yo ~ x + fc) + bf(yo2 ~ z + fc), family = sratio(), data = d)))
for (tag in c("WITH cs", "WITHOUT cs")) {
  f <- if (tag == "WITH cs") mv else mvn
  cat("\n[", tag, "]\n", sep = "")
  go("fitted() in sample", dim(fitted(f)))
  go("frm_linpred(resp='yo', newdata)",
     dim(frm_linpred(f, newdata = nd, resp = "yo", type = "response")))
  go("predict(resp='yo', newdata)",
     dim(predict(f, newdata = nd, resp = "yo", ndraws = 50)))
  go("simulate(nsim=2)", dim(as.data.frame(simulate(f, nsim = 2))))
  go("conditional_effects(effects='fc')",
     nrow(as.data.frame(conditional_effects(f, effects = "fc",
                                            resp = "yo")[[1L]])))
  go("summary() fixef rows", nrow(fixef(f)))
}
cat("\nDONE ", TAG, "\n")
