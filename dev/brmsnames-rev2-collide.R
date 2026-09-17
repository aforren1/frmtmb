## Reviewer recheck, MAJOR 2: collision classes the round-1 fix may not
## reach. (1) An interaction group whose levels carry brms's join
## character: `(1 | gi:hi)` with levels 1_2:3 and 1:2_3, which brms's
## combine_groups() pastes to the same level, so brms fits ONE group
## there. (2) Factor levels that brms's rename() maps to one name.
## (3) Coefficient names that rename() maps to one name.
##   Rscript dev/brmsnames-rev2-collide.R base|lane
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
show <- function(label, v) {
  if (is.character(v) && length(v) == 1L) {
    cat(sprintf("  %-34s %s\n", label, substr(gsub("\n", " ", v), 1, 160)))
  } else {
    cat(sprintf("  %-34s ", label)); print(v)
  }
}
q(library(frmtmb)); q(library(frmtmb.sample))
source("dev/brmsnames-rev2-data.R")
d <- rev2_data()
cat("arm", arm, "\n")

cat("\n== C1 (1 | gi:hi), levels 1_2:3 and 1:2_3 ==\n")
fit <- try1(q(frm(bf(y ~ x + (1 | gi:hi)), family = gaussian(), data = d)))
if (!is.character(fit)) {
  show("ranef(fit) rownames", try1(rownames(ranef(fit)[[1]])))
  show("coef(fit) rownames", try1(rownames(coef(fit)[[1]])))
  set.seed(3)
  ds <- try1(q(frm_sample(fit, chains = 1, iter = 100, refresh = 0,
                          seed = 3)))
  if (!is.character(ds)) {
    show("draws colnames", colnames(ds$draws))
    show("as_draws_df", try1(dim(as_draws_df(ds))))
    show("summary", try1(dim(summary(ds))))
    show("posterior_summary", try1(dim(posterior_summary(ds))))
    show("fixef(ds)", try1(dim(fixef(ds))))
    r <- try1(ranef(ds))
    show("ranef(ds) level names", if (is.list(r)) dimnames(r[[1]])[[1]] else r)
    show("coef(ds)", try1(dimnames(coef(ds)[[1]])[[1]]))
    show("VarCorr(ds)", try1(names(VarCorr(ds))))
    show("hypothesis on r_ name", try1(q(hypothesis(ds,
      "r_gi:hi[1_2_3,Intercept] = 0", class = NULL))$hypothesis$Estimate))
    show("column means of the two", try1(colMeans(ds$draws[,
      grep("1_2_3", colnames(ds$draws))])))
    show("log_lik dim", try1(dim(log_lik(ds))))
  }
}

cat("\n== C2 factor levels 'a b' and 'ab' (both fab in brms) ==\n")
# a reference level that sorts first, so both colliding levels get a
# column
d$fab3 <- factor(ifelse(seq_len(nrow(d)) %% 3 == 0, "0",
                        as.character(d$fab)), levels = c("0", "a b", "ab"))
fit <- try1(q(frm(bf(y ~ fab3), family = gaussian(), data = d)))
show("frm", if (is.character(fit)) fit else try1(variables(fit)))

cat("\n== C3 I(x^2) beside a column named IxE2 ==\n")
fit <- try1(q(frm(bf(y ~ I(x^2) + IxE2), family = gaussian(), data = d)))
show("frm", if (is.character(fit)) fit else try1(variables(fit)))

cat("\n== C4 group levels 'lvl 1' and 'lvl.1' in one factor ==\n")
gd <- as.character(d$g)
i1 <- which(gd == "lvl 1")
gd[i1[c(TRUE, FALSE)]] <- "lvl.1"
d$gd <- factor(gd)
fit <- try1(q(frm(bf(y ~ x + (1 | gd)), family = gaussian(), data = d)))
if (!is.character(fit)) {
  lab <- if (exists("brms_par_labels", asNamespace("frmtmb"))) frmtmb:::brms_par_labels(fit) else character(0)
  show("draws labels duplicated", lab[duplicated(lab)])
  set.seed(3)
  ds <- try1(q(frm_sample(fit, chains = 1, iter = 100, refresh = 0,
                          seed = 3)))
  if (!is.character(ds)) {
    show("as_draws_df", try1(dim(as_draws_df(ds))))
    r <- try1(ranef(ds))
    show("ranef(ds) level names",
         if (is.list(r)) dimnames(r[[1]])[[1]] else r)
  }
} else show("frm", fit)

cat("\n== C5 mv responses y_a and ya (both ya in brms) ==\n")
fit <- try1(q(frm(mvbf(bf(y_a ~ x), bf(ya ~ x), rescor = FALSE),
                  family = gaussian(), data = d)))
show("frm", if (is.character(fit)) fit else try1(variables(fit)))

cat("\n== C6 smooth by-factor levels 'u' in s(x, by = f2) beside s(x) ==\n")
fit <- try1(q(frm(bf(y ~ s(x) + s(z, by = f2)), family = gaussian(),
                  data = d)))
show("frm", if (is.character(fit)) fit else try1(variables(fit)))
