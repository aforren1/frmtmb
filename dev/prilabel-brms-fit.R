# Lane prilabel: one tiny brms 2.23.0 fit to record where the nested
# and multi-membership group labels surface: variables(), ranef(),
# coef(), VarCorr(), summary(), prior_summary(), and predict() with
# re_formula and new levels. Data: dev/dblbar-data.R, seed 20261013;
# sampling seed 1. The draws are not used, only the names.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(brms))
dd <- dblbar_data()
sink(file.path(wt, "dev/prilabel-log/brms-fit.txt"), split = TRUE)
cat("brms", format(packageVersion("brms")), "\n")
fit <- brm(bf(y ~ x + (1 + x | g1/g2) + (1 | mm(m1, m2)), sigma ~ (1 | g)),
           data = dd, chains = 1, iter = 300, warmup = 150, seed = 1,
           refresh = 0,
           prior = set_prior("normal(0, 0.5)", class = "sd",
                             group = "g1:g2", coef = "x"))
saveRDS(fit, file.path(wt, "dev/prilabel-log/brms-fit.rds"))
v <- variables(fit)
cat("\n== variables (sd_, cor_, first r_ of each group)\n")
print(grep("^sd_|^cor_", v, value = TRUE))
r <- grep("^r_", v, value = TRUE)
print(r[!duplicated(sub("\\[.*", "", r))])
print(head(grep("^r_g1:g2", r, value = TRUE), 3))
print(head(grep("^r_mm", r, value = TRUE), 3))
cat("\n== names(ranef), dimnames heads\n")
re <- ranef(fit)
print(names(re))
for (nm in names(re)) print(head(dimnames(re[[nm]])[[1]], 3))
print(dimnames(re[["g1:g2"]])[[3]])
cat("\n== names(coef)\n")
print(names(coef(fit)))
cat("\n== names(VarCorr)\n")
print(names(VarCorr(fit)))
cat("\n== summary random part\n")
s <- summary(fit)
print(names(s$random))
for (nm in names(s$random)) print(rownames(s$random[[nm]]))
cat("\n== prior_summary sd rows\n")
ps <- as.data.frame(prior_summary(fit))
print(ps[ps$class %in% c("sd", "cor"),
         c("prior", "class", "coef", "group", "dpar")], row.names = FALSE)
cat("\n== ngrps\n")
print(ngrps(fit))
cat("\n== predict, re_formula naming g1:g2 only\n")
nd <- dd[1:3, ]
p1 <- fitted(fit, newdata = nd, re_formula = ~ (1 + x | g1:g2))
print(dim(p1))
p2 <- tryCatch(fitted(fit, newdata = nd, re_formula = ~ (1 + x | g1/g2)),
               error = function(e) paste("ERROR:", conditionMessage(e)))
print(if (is.character(p2)) p2 else dim(p2))
p3 <- tryCatch(fitted(fit, newdata = nd, re_formula = ~ (1 + x | g2:g1)),
               error = function(e) paste("ERROR:", conditionMessage(e)))
print(if (is.character(p3)) p3 else dim(p3))
cat("\n== newdata with a new level of g2 within g1\n")
nd2 <- dd[1:2, ]
nd2$g2 <- factor("qNEW")
p4 <- tryCatch(fitted(fit, newdata = nd2, allow_new_levels = TRUE),
               error = function(e) paste("ERROR:", conditionMessage(e)))
print(if (is.character(p4)) p4 else dim(p4))
p5 <- tryCatch(fitted(fit, newdata = nd2),
               error = function(e) paste("ERROR:", conditionMessage(e)))
print(p5)
sink()
