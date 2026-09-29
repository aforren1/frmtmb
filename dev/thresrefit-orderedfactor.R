## The response as an ORDERED FACTOR, not integer codes: the pin makes
## the refit's family carry a threshold count, and frame.R refuses a
## count that asks for more categories than the response has levels.
## Does a leave-one-out subset drop the top level?
arm <- Sys.getenv("FRMTMB_LIB", "lane")
lib <- if (identical(arm, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("ARM", arm, "frmtmb", as.character(packageVersion("frmtmb")), "\n")
say <- function(...) cat(..., "\n", sep = "")

set.seed(501)
n <- 50
x <- rnorm(n)
cp <- cbind(plogis(-0.7 - 0.5 * x), plogis(0.6 - 0.5 * x),
            plogis(2.3 - 0.5 * x))
y <- 1L + rowSums(runif(n) > cp)
top <- which(y == 4L)
y[top[-1L]] <- 3L
dd <- data.frame(x = x, y = factor(y, levels = 1:4, ordered = TRUE))
itop <- which(as.integer(dd$y) == 4L)
say("table(y) = ", paste(table(dd$y), collapse = "/"),
    "; the single top row is ", itop)
say("levels kept by the subset: ",
    paste(levels(dd$y[-itop]), collapse = ","))

fit <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = dd))
say("fit tau_raw length = ", length(fit$estimates$tau_raw))
inf <- tryCatch(suppressWarnings(influence(fit, force = TRUE)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(inf)) say("influence -> ", inf) else {
  say("colnames = ", paste(colnames(inf$fixed), collapse = ","))
  say("NA cells = ", sum(is.na(inf$fixed)))
  say("row ", itop, " = ",
      paste(format(inf$fixed[itop, ], digits = 10), collapse = " "))
  cd <- cooks.distance(inf)
  say("cooks.distance NA count = ", sum(is.na(cd)),
      "; at row ", itop, " = ", format(cd[itop], digits = 8))
}

## and with the unused level DROPPED from the stored data, which is what
## a user who called droplevels() would hand influence(data = )
dd2 <- dd
dd2$y <- droplevels(dd2$y[drop = FALSE])
say("dropped-level variant: levels = ", paste(levels(dd2$y), collapse = ","))
