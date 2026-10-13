# Punch round 1, m7: what a || term with mo() fits as.
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
n <- 300
dd <- data.frame(x = sample(0:3, n, TRUE),
                 f = factor(sample(c("a", "b", "c"), n, TRUE)),
                 g = factor(sample(paste0("g", 1:20), n, TRUE)))
dd$y <- 1 + dd$x * (1 + (dd$f == "b")) + rnorm(20)[dd$g] + rnorm(n)
print(frmtmb:::expand_double_verts(y ~ mo(x) * f + (mo(x) * f || g)))
fit <- suppressWarnings(frm(bf(y ~ mo(x) * f + (mo(x) * f || g)),
                            family = gaussian(), data = dd))
print(vapply(fit$frame$re_blocks, function(b) b$covstruct, ""))
print(grep("^(sd|cor)_g", variables(fit), value = TRUE))
print(VarCorr(fit)$g)
f2 <- suppressWarnings(frm(bf(y ~ f + (f || g)), family = gaussian(),
                           data = dd))
print(vapply(f2$frame$re_blocks, function(b) b$covstruct, ""))
print(grep("^(sd|cor)_g", variables(f2), value = TRUE))
