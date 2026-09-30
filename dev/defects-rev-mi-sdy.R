# Reviewer of lane defects, recheck: mi(sd = ) measured rows read the
# measurement; subset() with mi() in a multivariate formula.
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("dev/defects-rev-mi-data.R")
d <- mi_data()
f <- frm(ymeas | mi(sdy) ~ x, data = d)
r <- residuals(f)
ok <- !is.na(d$ymeas)
cat("measured rows: residual == ymeas - fitted:",
    isTRUE(all.equal(unname(r[ok, "Estimate"]), d$ymeas[ok] - unname(fitted(f)[ok, "Estimate"]))), "\n")
cat("Est.Error at unmeasured rows (brms: NA):", r[!ok, "Est.Error"][1:3], "\n")
d$s <- d$x > -1
res <- tryCatch(frm(bf(y | subset(s) ~ x) + bf(xmi | mi() ~ z) + set_rescor(FALSE), data = d),
                error = function(e) conditionMessage(e))
cat("subset() + mi():", if (is.character(res)) substr(res, 1, 200) else "fitted", "\n")
