## Why does the pinned leave-one-out refit fail on an ORDERED FACTOR
## response? Call the seam directly so the error is visible.
.libPaths(c("C:/Users/adf44/source/r/wt-thresrefit-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
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
fit <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = dd))
pin <- frmtmb:::thres_pin_of_fit(fit)
say("pin: grouped = ", pin$y$grouped, " nthres = ", pin$y$nthres)
d0 <- fit$frame[["data_frame"]]
say("stored data_frame: class(y) = ", paste(class(d0$y), collapse = "/"),
    " levels = ", paste(levels(d0$y), collapse = ","))
sub <- d0[-itop, , drop = FALSE]
say("subset levels = ", paste(levels(sub$y), collapse = ","),
    " observed = ", paste(sort(unique(as.integer(sub$y))), collapse = ","))
r <- tryCatch(frmtmb:::assemble_frame(fit$spec, sub, thres_pin = pin),
              error = function(e) paste("ERROR:", conditionMessage(e)))
say("assemble_frame with pin -> ",
    if (is.character(r)) r else
      paste("tau_raw length", length(r$par_template$tau_raw)))
r2 <- tryCatch(frmtmb:::assemble_frame(fit$spec, sub),
               error = function(e) paste("ERROR:", conditionMessage(e)))
say("assemble_frame without pin -> ",
    if (is.character(r2)) r2 else
      paste("tau_raw length", length(r2$par_template$tau_raw)))
## what y_levels the frame records
if (!is.character(r2)) {
  say("y_levels without pin = ",
      paste(r2$y_levels$y, collapse = ","))
}
