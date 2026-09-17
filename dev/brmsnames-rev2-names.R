## Reviewer recheck, MAJOR 1 and 2 on hostile names: the lane's names
## (variables(fit), the draws labels) against brms's own fitted names
## (dev/brmsnames-rev2-brms.R), and hypothesis() on those names against a
## hand computation that uses neither parser.
##   Rscript dev/brmsnames-rev2-names.R base|lane
## Data seed 52 (dev/brmsnames-rev2-data.R), sampler seed 3.
arm <- commandArgs(trailingOnly = TRUE)[1L]
libs <- list(
  base = c("C:/Users/adf44/source/r/rellib-r3"),
  lane = c("C:/Users/adf44/source/r/brmsnames-lib",
           "C:/Users/adf44/source/r/rellib-r3"))
.libPaths(c(libs[[arm]], "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
try1 <- function(e) {
  tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
}
q(library(frmtmb)); q(library(frmtmb.sample))
source("dev/brmsnames-rev2-data.R")
d <- rev2_data()
ms <- rev2_models()
lane <- arm == "lane"
fam_of <- function(nm) {
  ns <- asNamespace("frmtmb")
  if (exists(nm, ns)) get(nm, ns)() else get(nm, asNamespace("stats"))()
}
internal <- "^(theta_|thetaac|thetar|b\\[|miss|lp__|lprior|Intercept)"
for (k in names(ms)) {
  M <- ms[[k]]
  cat("\n==", k, "==\n")
  f <- rev2_formula(k, M, frmtmb::bf, frmtmb::mvbf)
  fit <- try1(q(frm(f, family = fam_of(M$fam), data = d)))
  if (is.character(fit)) {
    cat("frm:", substr(fit, 1, 200), "\n")
    next
  }
  vf <- try1(variables(fit))
  lab <- try1(frmtmb:::brms_par_labels(fit))
  bf_path <- sprintf("dev/stan-cache/brmsnames-rev2-brms-%s.rds", k)
  bv <- if (file.exists(bf_path)) names(readRDS(bf_path)$fit@sim$samples[[1]]) else NULL
  cat("variables(fit):", vf, "\n")
  if (!is.character(lab) || length(lab) > 1L) {
    dup <- lab[duplicated(lab)]
    cat("draws labels duplicated:", if (length(dup)) dup else "none", "\n")
  } else cat("draws labels:", lab, "\n")
  if (!is.null(bv)) {
    lb <- if (length(lab) > 1L) lab else vf
    cat("lane labels not in brms (internal names dropped):",
        setdiff(grep(internal, lb, value = TRUE, invert = TRUE), bv), "\n")
    bn <- grep("^(lprior|lp__|Intercept|s_|zgp_|z_|simo_|L_|zb)", bv,
               value = TRUE, invert = TRUE)
    cat("brms names not in lane labels:", setdiff(bn, lb), "\n")
  }
  if (k == "B8") {
    bk <- fit$frame$re_blocks[[1]]
    cat("lane levels of gi:hi:", bk$levels, "\n")
  }
}

## hypothesis() attacks on B1, on the fit and on draws
cat("\n== hypothesis attacks, B1 ==\n")
fit <- q(frm(bf(ms$B1$f), family = gaussian(), data = d))
fe <- fixef(fit)
fe <- if (is.list(fe)) unlist(fe) else fe
names(fe) <- sub("^mu[.]", "", names(fe))
vcm <- if (lane) varcorr_matrices(fit) else VarCorr(fit)
sdx <- sqrt(diag(vcm[[1]]))[["x"]]
pick <- function(nm) fe[[nm]]
hand <- c(
  `x:z:fe:f > 0` = pick("x:z:fe:f"),
  `x:fcMd - IxE2 = 0` = pick("x:fc-d") - pick("I(x^2)"),
  `polyz2rawEQTRUE2 = 0` = pick("poly(z, 2, raw = TRUE)2"),
  `fe:f + x:fe:f = 0` = pick("fe:f") + pick("x:fe:f"),
  `b_x:z - sd_g__x > 0` = pick("x:z") - sdx
)
hs <- names(hand)
for (i in seq_along(hs)) {
  h <- hs[i]
  cls <- if (grepl("sd_g", h)) NULL else "b"
  r <- try1(if (lane) {
    q(hypothesis(fit, h, class = cls))$hypothesis$Estimate
  } else q(hypothesis(fit, h))$estimate)
  cat(sprintf("fit  %-24s hand %.6f  hypothesis %s\n", h, hand[[i]],
              paste(r, collapse = " ")))
}
if (lane) {
  r <- q(hypothesis(fit, hs[1:4]))$hypothesis$Estimate
  cat("fit  four at once:", sprintf("%.6f", r), " hand:",
      sprintf("%.6f", hand[1:4]), "\n")
  set.seed(3)
  ds <- q(frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 3))
  a <- as_draws_df(ds)
  sdraw <- VarCorr(ds, summary = FALSE)$g$sd[, "x"]
  hd <- c(mean(a[["b_x:z:fe:f"]]),
          mean(a[["b_x:fcMd"]] - a[["b_IxE2"]]),
          mean(a[["b_polyz2rawEQTRUE2"]]),
          mean(a[["b_fe:f"]] + a[["b_x:fe:f"]]),
          mean(a[["b_x:z"]] - sdraw))
  for (i in seq_along(hs)) {
    cls <- if (grepl("sd_g", hs[i])) NULL else "b"
    r <- try1(q(hypothesis(ds, hs[i], class = cls))$hypothesis$Estimate)
    cat(sprintf("draws %-24s hand %.6f  hypothesis %s\n", hs[i], hd[i],
                paste(r, collapse = " ")))
  }
  co <- coef(ds, summary = FALSE)$g
  lv <- dimnames(co)[[2]]
  hc <- vapply(lv, function(l) mean(co[, l, "Intercept"] + co[, l, "x:fe:f"]),
               1)
  r <- q(hypothesis(ds, "Intercept + x:fe:f > 0", scope = "coef",
                    group = "g"))$hypothesis
  cat("draws scope coef group g, Intercept + x:fe:f > 0: max |hand - est|",
      max(abs(hc - r$Estimate[match(lv, r$Group)])), " levels",
      length(lv), " groups in result", nlevels(r$Group), "\n")
  r2 <- try1(q(hypothesis(ds, "r_g[lvl.1,x] - r_g[lvl.2,x] = 0",
                          class = NULL))$hypothesis$Estimate)
  cat("draws r_g[lvl.1,x] - r_g[lvl.2,x]: hand",
      mean(a[["r_g[lvl.1,x]"]] - a[["r_g[lvl.2,x]"]]), " hypothesis", r2,
      "\n")
}
