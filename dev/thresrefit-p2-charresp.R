## Nit 3, verified before it is filed: a character-coded ordinal response
## whose labels are not numeric text. Both arms.
arm <- Sys.getenv("FRMTMB_LIB", "lane")
lib <- if (identical(arm, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("ARM", arm, "frmtmb", as.character(packageVersion("frmtmb")), "\n")
say <- function(...) cat(..., "\n", sep = "")
msg <- function(expr) {
  tryCatch({ expr; "NO ERROR" }, error = function(e) conditionMessage(e))
}

set.seed(2501)
n <- 60
x <- rnorm(n)
cp <- cbind(plogis(-0.7 - 0.5 * x), plogis(0.6 - 0.5 * x),
            plogis(2.3 - 0.5 * x))
y <- 1L + rowSums(runif(n) > cp)

## numeric text works, as the docs say
d_num <- data.frame(x = x, y = as.character(y))
say("character, numeric labels ('1'..'4') -> ",
    msg(suppressWarnings(frm(bf(y ~ x), family = cumulative(),
                             data = d_num))))

## non-numeric labels do not
lab <- c("none", "mild", "moderate", "severe")
d_chr <- data.frame(x = x, y = lab[y])
say("character, non-numeric labels -> ",
    msg(suppressWarnings(frm(bf(y ~ x), family = cumulative(),
                             data = d_chr))))

## the ordered factor with the SAME labels is the supported spelling
d_ord <- data.frame(x = x, y = factor(lab[y], levels = lab, ordered = TRUE))
f <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = d_ord))
say("ordered factor, the same labels -> fits, n_tau = ",
    length(f$estimates$tau_raw))
