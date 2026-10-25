# Reviewer: frmtmb's brmsfamily("cox", bhaz =) and the compat rows of bhaz().
.libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb.survival))
r <- tryCatch(brmsfamily("cox", bhaz = list(df = 6)), error = function(e) conditionMessage(e))
cat("frmtmb brmsfamily('cox', bhaz =):", if (is.character(r)) r else "returned a family", "\n")
cc <- frm_compat("bhaz()")
cat("bhaz() compat rows:", nrow(cc), " status table:\n"); print(table(cc$status))
print(cc[cc$feature_b %in% c("cox", "royston_parmar", "gaussian"), c("feature_b", "status")])
print(cc[cc$status %in% c("works", "conditional", "untested"), c("feature_b", "kind_b", "status")], row.names = FALSE)
set.seed(1); d <- data.frame(time = rexp(50), cens = 0, x = rnorm(50))
r <- tryCatch(frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = royston_parmar(), data = d), error = function(e) conditionMessage(e))
cat("royston_parmar + bhaz:", if (is.character(r)) r else "FIT", "\n")
